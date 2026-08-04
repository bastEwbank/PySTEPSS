import pandapower as pp
import numpy as np
import plotly.graph_objects as go
from pyproj import Transformer
from plotly.colors import sample_colorscale
from pathlib import Path
import warnings

def normalize_colorscale(cmap):
    """
    Normalize the user-provided 'cmap' argument into a format that Plotly can use
    for both:
      - `sample_colorscale` (map a value -> a color)
      - the `marker.colorscale` attribute of a trace (for the colorbar)

    Parameters
    ----------
    cmap : str, list, or tuple
        - str: Name of a known Plotly colorscale (e.g., 'Jet', 'RdBu', 'Viridis').
        - list/tuple of colors, e.g., ['blue', 'white', 'red'] -> evenly distributed.
        - list/tuple of [fraction, color] pairs already in Plotly format,
          e.g., [[0, 'blue'], [0.5, 'white'], [1, 'red']] -> used as-is (custom cmap).

    Returns
    -------
    str or list
        Normalized colorscale compatible with Plotly.

    Raises
    ------
    ValueError
        If `cmap` is an empty list or an unsupported type.
    """
    if isinstance(cmap, str):
        return cmap

    if isinstance(cmap, (list, tuple)):
        if len(cmap) == 0:
            raise ValueError("cmap cannot be an empty list")

        first = cmap[0]
        # Already in [fraction, color] format?
        if isinstance(first, (list, tuple)) and len(first) == 2 and \
           isinstance(first[0], (int, float)):
            return [list(pair) for pair in cmap]

        # Otherwise, simple list of colors -> distribute evenly over [0,1]
        n = len(cmap)
        if n == 1:
            return [[0, cmap[0]], [1, cmap[0]]]
        return [[i / (n - 1), color] for i, color in enumerate(cmap)]

    raise ValueError(
        "cmap must be either a Plotly colorscale name (str), "
        "a list of colors, or a list of [fraction, color] pairs."
    )

def get_color_for_value(value, vmin, vmax, colorscale):
    """
    Return the color (str) corresponding to 'value' according to the given
    colorscale, normalized between vmin and vmax (clamped to [0,1]).

    Parameters
    ----------
    value : float or None
        The value to map to a color.
    vmin : float
        Minimum value of the range.
    vmax : float
        Maximum value of the range.
    colorscale : str or list
        The colorscale to use (normalized via `normalize_colorscale`).

    Returns
    -------
    str
        The color corresponding to the normalized value.
    """
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return 'grey'

    if vmax == vmin:
        norm = 0.5
    else:
        norm = (value - vmin) / (vmax - vmin)
    norm = min(max(norm, 0.0), 1.0)

    return sample_colorscale(colorscale, norm)[0]

def compute_bounding_box(coordinates: np.array):
    """
    Compute the width and height of the bounding box for a set of 2D coordinates.

    Parameters
    ----------
    coordinates : np.ndarray
        Array of shape (N, 2) containing the x-y or lon-lat coordinates.

    Returns
    -------
    tuple
        (width, height) of the bounding box.
    """
    x = coordinates[:, 0]
    y = coordinates[:, 1]

    min_x, max_x = min(x), max(x)
    min_y, max_y = min(y), max(y)

    return max_x - min_x, max_y - min_y

def compute_bus_marker_size(p_mw_res, p_mw_max, size_min=6, size_max=30):
    """
    Compute the marker size for a bus based on the magnitude of its net active
    power (production or consumption). The size is capped at `size_max` once
    |p_mw_res| reaches `p_mw_max`. At p_mw_res = 0, the size is `size_min`.

    Parameters
    ----------
    p_mw_res : float
        Net active power at the bus (MW). Positive for consumption,
        negative for production.
    p_mw_max : float or None
        Maximum absolute net active power for scaling. If None or <= 0,
        `size_min` is returned.
    size_min : float, optional
        Minimum marker size (default: 6).
    size_max : float, optional
        Maximum marker size (default: 30).

    Returns
    -------
    float
        The computed marker size.
    """
    if p_mw_max is None or p_mw_max <= 0:
        return size_min

    ratio = min(abs(p_mw_res) / p_mw_max, 1.0)
    return size_min + ratio * (size_max - size_min)

def compute_bus_marker_symbol(p_mw_res, eps=1e-6):
    """
    Return the marker symbol for a bus depending on the sign of `p_mw_res`.
    In pandapower convention:
      - p_mw_res > 0 -> net consumption (hollow circle: 'circle-open').
      - p_mw_res <= 0 -> net production or balanced (filled circle: 'circle').

    Parameters
    ----------
    p_mw_res : float
        Net active power at the bus (MW).
    eps : float, optional
        Threshold for considering a value as positive (default: 1e-6).

    Returns
    -------
    str
        The marker symbol ('circle-open' or 'circle').
    """
    if p_mw_res > eps:
        return 'circle-open'
    return 'circle'

def append_arrow_trace(p_from: np.array, p_to: np.array, arrow_traces, on_map=False,
                       arrow_color='blue', arrow_scaling=None, direction='from_to'):
    """
    Append an arrow trace to `arrow_traces` to represent power flow direction
    on a line or transformer. The arrow is placed at the end of the edge.

    Parameters
    ----------
    p_from : np.ndarray
        Starting point coordinates (x, y) or (lon, lat).
    p_to : np.ndarray
        Ending point coordinates (x, y) or (lon, lat).
    arrow_traces : list
        List to which the arrow trace will be appended.
    on_map : bool, optional
        If True, coordinates are treated as lon-lat and transformed to Web
        Mercator (EPSG:3857) for arrow geometry calculations (default: False).
    arrow_color : str, optional
        Color of the arrow (default: 'blue').
    arrow_scaling : float or None, optional
        Scaling factor for arrow size. If None, a default scaling is used.
    direction : str, optional
        Direction of the arrow: 'from_to' or 'to_from' (default: 'from_to').
    """
    # Compute direction vector
    if direction == 'to_from':
        v = p_from - p_to
    else:  # 'from_to' + other
        v = p_to - p_from

    v_length = np.linalg.norm(v)
    if v_length == 0:
        return  # Avoid division by zero
    w = v / v_length

    # Default or custom arrow length
    if arrow_scaling is None:
        h = v_length * 0.05  # Default: 5% of edge length
    else:
        h = (0.02 * arrow_scaling + 0.001) * (1 - np.exp(-12 * v_length / arrow_scaling))

    center_edge = (p_from + p_to) / 2 + w * h * 0.5  # Shifted center
    P = center_edge - h * w  # Bottom of the arrow

    if on_map:
        # Transform to Web Mercator for accurate geometry
        transformer = Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True)
        x1, y1 = transformer.transform(P[0], P[1])
        x2, y2 = transformer.transform(center_edge[0], center_edge[1])
        width = h * 1e5  # Adjust scaling for projected coordinates
    else:
        x1, y1 = P[0], P[1]
        x2, y2 = center_edge[0], center_edge[1]
        width = h / 2

    P_prime = np.array([x1, y1])
    center_edge_prime = np.array([x2, y2])
    v_prime = center_edge_prime - P_prime
    w_prime = v_prime / np.linalg.norm(v_prime)
    u_prime = np.array([-w_prime[1], w_prime[0]])  # Orthogonal vector

    S_prime = P_prime - width * u_prime
    T_prime = P_prime + width * u_prime

    if on_map:
        # Transform back to lon-lat
        transformer_reverse = Transformer.from_crs("EPSG:3857", "EPSG:4326", always_xy=True)
        S = transformer_reverse.transform(S_prime[0], S_prime[1])
        T = transformer_reverse.transform(T_prime[0], T_prime[1])
        go_func = go.Scattermap
        selected_params = ['lon', 'lat']
    else:
        S = S_prime
        T = T_prime
        go_func = go.Scatter
        selected_params = ['x', 'y']

    args_pos = {key: [S[i], P[i], T[i], center_edge[i], S[i]]
                for i, key in enumerate(selected_params)}

    trace = go_func(
        mode='lines',
        fill='toself',
        fillcolor=arrow_color,
        line_color=arrow_color,
        opacity=1,
        hoverinfo='skip',
        **args_pos
    )
    arrow_traces.append(trace)

def compute_box_center(coordinates: np.array):
    """
    Compute the center of the bounding box for a set of 2D coordinates.

    Parameters
    ----------
    coordinates : np.ndarray
        Array of shape (N, 2) containing the x-y or lon-lat coordinates.

    Returns
    -------
    tuple
        (center_x, center_y) of the bounding box.
    """
    x = coordinates[:, 0]
    y = coordinates[:, 1]
    min_x, max_x = min(x), max(x)
    min_y, max_y = min(y), max(y)
    center_x = (min_x + max_x) / 2
    center_y = (min_y + max_y) / 2
    return center_x, center_y

def pf_res_plotly_custom(
    net: pp.pandapowerNet,
    on_map: bool = False,
    filename: str = None,
    map_style: str = 'basic',
    climits_volt: tuple = (0.9, 1.1),
    climits_lines: tuple = (0, 100),
    line_width: float = 2,
    bus_size: float = 10,
    bus_size_min: float = None,
    bus_size_max: float = None,
    bus_marker_line_width: float = 3,
    cpos_volt: float = 1.0,
    cpos_lines: float = 1.1,
    cmap_volt: str = 'Jet',
    cmap_lines: str = 'Jet',
) -> go.Figure:
    """
    Customize the plot generated by `pandapower.plotting.plotly.pf_res_plotly`
    by directly manipulating the `data` tuple of the `plotly.go.Figure`.

    This function enhances the default power flow result plot with:
      - Customizable colormaps for voltage and line loading.
      - Bus marker sizes scaled by net active power.
      - Bus marker symbols indicating production (filled) or consumption (hollow).
      - Arrows on lines and transformers indicating power flow direction.
      - Enhanced hover information for buses, lines, and transformers.
      - Support for both Cartesian and geographic (map) plots.

    Parameters
    ----------
    net : pp.pandapowerNet
        The pandapower network object.
    on_map : bool, optional
        If True, the plot is rendered on a map (lon-lat coordinates).
        If False, a Cartesian plot is used (default: False).
    filename : str, optional
        Output filename for saving the plot. Supported formats: .html, .png,
        .jpeg, .webp, .svg, .pdf (requires Kaleido for non-HTML formats).
        Default: 'temp_plot.html'.
    map_style : str, optional
        Style of the map if `on_map` is True (default: 'basic').
    climits_volt : tuple, optional
        Color limits for voltage magnitude (pu). Default: (0.9, 1.1).
    climits_lines : tuple, optional
        Color limits for line loading percent. Default: (0, 100).
    line_width : float, optional
        Width of the lines in the plot (default: 2).
    bus_size : float, optional
        Base size for bus markers (default: 10).
    bus_size_min : float or None, optional
        Minimum bus marker size. If None, `bus_size` is used (default: None).
    bus_size_max : float or None, optional
        Maximum bus marker size. If None, `bus_size * 3` is used (default: None).
    bus_marker_line_width : float, optional
        Width of the bus marker outline (default: 3).
    cpos_volt : float, optional
        Center position for voltage colormap (default: 1.0).
    cpos_lines : float, optional
        Center position for line loading colormap (default: 1.1).
    cmap_volt : str, optional
        Colormap for voltage magnitude (default: 'Jet').
    cmap_lines : str, optional
        Colormap for line loading (default: 'Jet').

    Returns
    -------
    go.Figure
        The customized Plotly figure object.
    """
    # Generate the initial Plotly figure with pandapower
    fig = pp.plotting.plotly.pf_res_plotly(
        net,
        on_map=on_map,
        filename=filename,
        line_width=line_width,
        bus_size=bus_size,
        climits_volt=climits_volt,
        climits_load=climits_lines,
        cpos_volt=cpos_volt,
        cpos_load=cpos_lines,
        auto_open=False,
        cmap='Jet'
    )

    # Check if the figure is on a map or not
    if isinstance(fig.data[0], go.Scatter):
        on_map = False
    elif isinstance(fig.data[0], go.Scattermap):
        on_map = True
    else:
        raise ValueError(
            f"Plotly figure is neither of type go.Scatter or go.Scattermap. "
            f"Got: {type(fig.data[0])}"
        )

    # --- Constants and variables declaration ---
    if on_map:
        x_key = 'lon'
        y_key = 'lat'
    else:
        x_key = 'x'
        y_key = 'y'

    # Normalize the user-provided colormaps
    cmap_lines_norm = normalize_colorscale(cmap_lines)
    cmap_volt_norm = normalize_colorscale(cmap_volt)

    number_of_buses = net.bus.shape[0]
    number_of_lines = net.line.shape[0]
    number_of_trafos = net.trafo.shape[0]

    min_bus_idx = min(net.bus.index)
    max_bus_idx = max(net.bus.index)
    adj_shape = max_bus_idx - min_bus_idx + 1  # Shape of adjacency matrix
    adjacency_matrix = np.zeros((adj_shape, adj_shape), dtype=tuple)

    bus_coords = np.empty((number_of_buses, 2))  # x-y or lon-lat coordinates

    
    arrow_traces = []

    # Compute maximum net active power for bus marker scaling
    p_mw_load_max = net.res_bus['p_mw'].abs().max()
    p_mw_load_max = max(p_mw_load_max, net.load['p_mw'].abs().max())
    p_mw_load_max = max(p_mw_load_max, net.sgen['p_mw'].abs().max())
    p_mw_load_max = max(p_mw_load_max, net.gen['p_mw'].abs().max())

    # Resolve default min/max bus marker sizes
    if bus_size_min is None:
        bus_size_min = bus_size
    if bus_size_max is None:
        bus_size_max = bus_size * 3

    bus_marker_sizes = []    # Size per bus, scaled with |p_mw_res|
    bus_marker_symbols = []  # 'circle' (producing/balanced) vs 'circle-open' (consuming)
    bus_marker_colors = []

    # --- Update Bus Trace ---
    # Enhance hover with additional result information, record geo coords,
    # and change marker shape/size according to net.res_bus['p_mw']
    for hover_str, idx in zip(fig.data[-1]['text'], np.arange(number_of_buses)):
        bus_id = hover_str.split('<br />')[0]  # Get bus name from hover string
        bus_idx = net.bus.index[net.bus['name'] == bus_id][0]  # Get bus index

        
        
        bus_coords[idx, :] = np.array([
            fig.data[-1][x_key][idx],
            fig.data[-1][y_key][idx]
        ])

        p_mw_res = net.res_bus.at[bus_idx, 'p_mw']  # Net active power (MW)
        q_mvar_res = net.res_bus.at[bus_idx, 'q_mvar']  # Net reactive power (MVar)

        # Get load, gen, and sgen values for detailed hover
        sgen_idx = net.sgen.index[net.sgen['bus'] == bus_idx].tolist()
        gen_idx = net.gen.index[net.gen['bus'] == bus_idx].tolist()
        load_idx = net.load.index[net.load['bus'] == bus_idx].tolist()

        if len(load_idx) > 0:
            p_mw_res_load = net.res_load.loc[load_idx, 'p_mw'].sum()
            q_mvar_res_load = net.res_load.loc[load_idx, 'q_mvar'].sum()
            hover_str += (
                f'<br />Load P = {p_mw_res_load:.2f} MW'
                f'<br />Load Q = {q_mvar_res_load:.2f} MVar'
            )
        if len(gen_idx) > 0:
            p_mw_res_gen = net.res_gen.loc[gen_idx, 'p_mw'].sum()
            q_mvar_res_gen = net.res_gen.loc[gen_idx, 'q_mvar'].sum()
            hover_str += (
                f'<br />Gen P = {p_mw_res_gen:.2f} MW'
                f'<br />Gen Q = {q_mvar_res_gen:.2f} MVar'
            )
        if len(sgen_idx) > 0:
            p_mw_res_sgen = net.res_sgen.loc[sgen_idx, 'p_mw'].sum()
            q_mvar_res_sgen = net.res_sgen.loc[sgen_idx, 'q_mvar'].sum()
            hover_str += (
                f'<br />SGen P = {p_mw_res_sgen:.2f} MW'
                f'<br />SGen Q = {q_mvar_res_sgen:.2f} MVar'
            )

        # Add net power
        hover_str += f'<br />P = {p_mw_res:.2f} MW<br />Q = {q_mvar_res:.2f} MVar'

        # Add bus index from the net dataframe
        hover_str = f"{bus_id} <br />Index = {bus_idx} <br />" + \
                   hover_str.split('<br />', 1)[1]

        # Update the hover
        fig.data[-1]['text'][idx] = hover_str

        # Compute bus marker properties
        bus_computed_size =compute_bus_marker_size(p_mw_res, p_mw_load_max,
                                                   size_min=bus_size_min, 
                                                   size_max=bus_size_max)
        bus_marker_sizes.append(bus_computed_size)
        bus_marker_symbols.append(compute_bus_marker_symbol(p_mw_res))

        # Recolor bus based on voltage magnitude
        v_mag = net.res_bus.at[bus_idx, 'vm_pu']
        bus_color = get_color_for_value(
            v_mag, climits_volt[0], climits_volt[1], cmap_volt_norm
        )
        bus_marker_colors.append(bus_color)
        
        #workaround we assume only one slack node per simulation 
        if net.ext_grid.iloc[0]['bus'] == bus_idx:
            # Increase size of slack bus marker to be visible with new marker styles
            fig.data[-2]['marker']['size'] =bus_computed_size + 5  # TODO: Follow bus marker size with offset
            fig.data[-2]['hoverinfo'] = 'skip'
            

    # --- Recolor the voltage colorbar if present ---
    if getattr(fig.data[-1]['marker'], 'colorbar') is not None:
        volt_colorbar_trace = fig.data[-1]
        if 'marker' in volt_colorbar_trace:
            volt_colorbar_trace['marker']['colorscale'] = cmap_volt_norm
            if getattr(volt_colorbar_trace['marker'], 'cmin') is None:
                volt_colorbar_trace['marker']['cmin'] = climits_volt[0]
            if getattr(volt_colorbar_trace['marker'], 'cmax') is None:
                volt_colorbar_trace['marker']['cmax'] = climits_volt[1]

    # Apply the per-bus size/symbol/color arrays
    bus_marker = fig.data[-1]['marker']
    bus_marker['size'] = bus_marker_sizes
    bus_marker['opacity'] = 1
    bus_marker['color'] = bus_marker_colors

    if not on_map:
        bus_marker['symbol'] = bus_marker_symbols
        bus_marker['line'] = dict(
            color=bus_marker['color'],
            colorscale=bus_marker['colorscale'],
            cmin=bus_marker['cmin'],
            cmax=bus_marker['cmax'],
            width=bus_marker_line_width,
        )

    # --- Compute layout characteristics after recording bus coords ---
    hori, verti = compute_bounding_box(bus_coords)  # Bounding box dimensions
    diag_box = np.sqrt(hori**2 + verti**2)
    zoomlevel = 1 + 19 * np.exp(-0.3 * diag_box)  # Zoom level based on diagonal

    # --- Update Lines Trace ---
    # Enhance hover with additional result information, build adjacency matrix,
    # and recolor lines to match the requested color gradient
    trace_idx = 0
    for line_data in fig.data[0:number_of_lines]:
        hover_str = line_data['text']
        line_id = hover_str.split('<br />')[0]  # Get line name from hover string
        line_idx = net.line.index[net.line['name'] == line_id][0]  # Get line index

        bus_from_idx = net.line.at[line_idx, 'from_bus']
        bus_to_idx = net.line.at[line_idx, 'to_bus']

        adj_from_idx = bus_from_idx - min_bus_idx
        adj_to_idx = bus_to_idx - min_bus_idx

        if adjacency_matrix[adj_from_idx, adj_to_idx] == 0:
            adjacency_matrix[adj_from_idx, adj_to_idx] = (trace_idx,)
        else:
            adjacency_matrix[adj_from_idx, adj_to_idx] += (trace_idx,)

        # Add line length in km to hover string
        line_length_km = net.line.at[line_idx, 'length_km']
        hover_str += f'Length = {line_length_km:.2f} km<br />'

        p_from_mw = net.res_line.at[line_idx, 'p_from_mw']
        p_to_mw = net.res_line.at[line_idx, 'p_to_mw']
        q_from_mvar = net.res_line.at[line_idx, 'q_from_mvar']
        q_to_mvar = net.res_line.at[line_idx, 'q_to_mvar']

        p_losses_mw = p_from_mw + p_to_mw
        q_losses_mvar = q_from_mvar + q_to_mvar
        p_flow = abs(p_from_mw) - abs(p_to_mw)

        if p_flow < 0:
            hover_str += (
                f'P to = {p_to_mw:.2f} MW'
                f'<br />Q to = {q_to_mvar:.2f} MVar'
                f'<br />P losses = {p_losses_mw:.2f} MW'
                f'<br />Q losses = {q_losses_mvar:.2f} MVar'
                f'<br />P from = {p_from_mw:.2f} MW'
                f'<br />Q from = {q_from_mvar:.2f} MVar<br />'
            )
        else:
            hover_str += (
                f'P from = {p_from_mw:.2f} MW'
                f'<br />Q from = {q_from_mvar:.2f} MVar'
                f'<br />P losses = {p_losses_mw:.2f} MW'
                f'<br />Q losses = {q_losses_mvar:.2f} MVar'
                f'<br />P to = {p_to_mw:.2f} MW'
                f'<br />Q to = {q_to_mvar:.2f} MVar<br />'
            )

        # Add line index and from/to bus names
        hover_str = (
            f"{line_id} <br />Index = {line_idx} <br />"
            f"From = {net.bus.loc[bus_from_idx, 'name']} <br />"
            f"To = {net.bus.loc[bus_to_idx, 'name']} <br />" +
            hover_str.split('<br />', 1)[1]
        )
        line_data['text'] = hover_str
        if line_data['showlegend']:
            line_data['showlegend'] = False

        # Recolor the line based on loading_percent
        loading_percent = net.res_line.at[line_idx, 'loading_percent']
        line_color = get_color_for_value(
            loading_percent, climits_lines[0], climits_lines[1], cmap_lines_norm
        )
        line_data['line']['color'] = line_color

        trace_idx += 1

    start_trafo = number_of_lines + 1  # +1 for the group of line traces
    trace_idx += 1  # Same reason

    # --- Recolor the loading percent colorbar if present ---
    color_bar = 0
    if getattr(fig.data[number_of_lines]['marker'], 'colorbar') is not None:
        start_trafo += 1
        trace_idx += 1  # +1 for the colorbar trace between lines and trafos
        color_bar = 1
        line_colorbar_trace = fig.data[number_of_lines]
        line_colorbar_trace['hoverinfo'] = 'skip'
        if 'marker' in line_colorbar_trace:
            line_colorbar_trace['marker']['colorscale'] = cmap_lines_norm
            if getattr(line_colorbar_trace['marker'], 'cmin') is None:
                line_colorbar_trace['marker']['cmin'] = climits_lines[0]
            if getattr(line_colorbar_trace['marker'], 'cmax') is None:
                line_colorbar_trace['marker']['cmax'] = climits_lines[1]

    # --- Update Transformers Trace ---
    # Enhance hover with additional result information, build adjacency matrix
    # for double trafos, and recolor trafos to match the requested color gradient
    for trafo_data in fig.data[start_trafo: start_trafo + number_of_trafos]:
        hover_str = trafo_data['text']
        trafo_id = hover_str.split('<br />')[0]  # Get trafo name from hover string
        trafo_idx = net.trafo.index[net.trafo['name'] == trafo_id][0]  # Get trafo index

        bus_from_idx = net.trafo.at[trafo_idx, 'hv_bus']
        bus_to_idx = net.trafo.at[trafo_idx, 'lv_bus']

        adj_from_idx = bus_from_idx - min_bus_idx
        adj_to_idx = bus_to_idx - min_bus_idx

        if adjacency_matrix[adj_from_idx, adj_to_idx] == 0:
            adjacency_matrix[adj_from_idx, adj_to_idx] = (trace_idx,)
        else:
            adjacency_matrix[adj_from_idx, adj_to_idx] += (trace_idx,)

        p_hv_mw = net.res_trafo.at[trafo_idx, 'p_hv_mw']
        p_lv_mw = net.res_trafo.at[trafo_idx, 'p_lv_mw']
        q_hv_mvar = net.res_trafo.at[trafo_idx, 'q_hv_mvar']
        q_lv_mvar = net.res_trafo.at[trafo_idx, 'q_lv_mvar']

        p_losses_mw = p_hv_mw + p_lv_mw
        q_losses_mvar = q_hv_mvar + q_lv_mvar
        p_flow = abs(p_hv_mw) - abs(p_lv_mw)

        if p_flow < 0:
            hover_str += (
                f'P LV = {p_lv_mw:.2f} MW<br />Q LV = {q_lv_mvar:.2f} MVar'
                f'<br />P losses = {p_losses_mw:.2f} MW'
                f'<br />Q losses = {q_losses_mvar:.2f} MVar<br />'
            )
        else:
            hover_str += (
                f'P HV = {p_hv_mw:.2f} MW<br />Q HV = {q_hv_mvar:.2f} MVar'
                f'<br />P losses = {p_losses_mw:.2f} MW'
                f'<br />Q losses = {q_losses_mvar:.2f} MVar<br />'
            )

        # Add trafo index from the net dataframe
        hover_str = f"{trafo_id} <br />Index = {trafo_idx} <br />" + \
                   hover_str.split('<br />', 1)[1]
        trafo_data['text'] = hover_str

        # Recolor the trafo based on loading_percent
        loading_percent = net.res_trafo.at[trafo_idx, 'loading_percent']
        trafo_color = get_color_for_value(
            loading_percent, climits_lines[0], climits_lines[1], cmap_lines_norm
        )
        trafo_data['line']['color'] = trafo_color

        trace_idx += 1

    # --- Handle multiple lines/transformers for the same bus pair ---
    # Add arrows to lines and trafos in the convention from-to / hv-lv
    where_idx = np.where(adjacency_matrix != 0)  # Non-zero entries
    for trace_tuple in adjacency_matrix[where_idx]:
        n_multi_lines = len(trace_tuple)
        edge_name = fig.data[trace_tuple[0]]['text'].split('<br />')[0].strip()

        # Get power flow direction to adapt arrow direction
        if edge_name in net.line['name'].values:
            edge_idx = net.line.index[net.line['name'] == edge_name][0]
            p_from = net.res_line.at[edge_idx, 'p_from_mw']
            p_to = net.res_line.at[edge_idx, 'p_to_mw']
            p_flow = abs(p_from) - abs(p_to)
        elif edge_name in net.trafo['name'].values:
            edge_idx = net.trafo.index[net.trafo['name'] == edge_name][0]
            p_hv = net.res_trafo.at[edge_idx, 'p_hv_mw']
            p_lv = net.res_trafo.at[edge_idx, 'p_lv_mw']
            p_losses = p_hv + p_lv
            p_flow = abs(p_hv) - abs(p_lv)
        else:
            p_flow = 0
            warnings.warn(f"Edge name {edge_name} not found in lines or transformers.")

        # Determine arrow direction based on power flow
        if p_flow < -10e-6:  # Small threshold to avoid floating point issues
            direction = 'to_from'
        elif p_flow > 10e-6:
            direction = 'from_to'
        else:  # Also when 0 losses
            if edge_name in net.line['name'].values:
                direction = 'from_to' if p_from >= 0 else 'to_from'
            elif edge_name in net.trafo['name'].values:
                direction = 'from_to' if p_hv >= 0 else 'to_from'

        if n_multi_lines > 1:  # Multiple traces for the same bus pair
            x_line = fig.data[trace_tuple[0]][x_key]
            y_line = fig.data[trace_tuple[0]][y_key]

            v_line = np.array([x_line[-1] - x_line[0], y_line[-1] - y_line[0]])
            v_length = np.linalg.norm(v_line)
            if v_length == 0:
                warnings.warn(
                    f"Zero length vector for trace {trace_tuple[0]} at bus pair {where_idx}. "
                    "Skipping arrow creation."
                )
                continue
            w_line = v_line / v_length
            u_line = np.array([-w_line[1], w_line[0]])  # Orthogonal vector

            p_center = np.array([x_line[1], y_line[1]])  # Center point of the line

            span_dist = 5e-2 * np.exp(-v_length)  # Inversely proportional spacing
            offset_dist = np.linspace(-span_dist, span_dist, n_multi_lines)

            for t_idx, o_dist in zip(trace_tuple, offset_dist):
                new_x = (
                    fig.data[t_idx][x_key][0],
                    p_center[0] + o_dist * u_line[0],
                    fig.data[t_idx][x_key][2]
                )
                new_y = (
                    fig.data[t_idx][y_key][0],
                    p_center[1] + o_dist * u_line[1],
                    fig.data[t_idx][y_key][2]
                )

                fig.data[t_idx][x_key] = new_x
                fig.data[t_idx][y_key] = new_y
                p_from = np.array([new_x[0], new_y[0]])
                p_to = np.array([new_x[1], new_y[1]])

                trace_color = fig.data[t_idx]['line']['color']
                append_arrow_trace(
                    p_from, p_to, arrow_traces,
                    arrow_color=trace_color,
                    on_map=on_map,
                    arrow_scaling=diag_box,
                    direction=direction
                )
        else:
            p_from = np.array([
                fig.data[trace_tuple[0]][x_key][0],
                fig.data[trace_tuple[0]][y_key][0]
            ])
            p_to = np.array([
                fig.data[trace_tuple[0]][x_key][2],
                fig.data[trace_tuple[0]][y_key][2]
            ])

            trace_color = fig.data[trace_tuple[0]]['line']['color']
            append_arrow_trace(
                p_from, p_to, arrow_traces,
                arrow_color=trace_color,
                on_map=on_map,
                arrow_scaling=diag_box,
                direction=direction
            )

    # --- Remove unused data traces of lines and trafos ---
    if number_of_trafos != 0 and number_of_lines != 0:
        fig.data = (
            fig.data[:number_of_lines + color_bar] +
            fig.data[start_trafo:start_trafo + number_of_trafos] +
            fig.data[start_trafo + number_of_trafos + 1:]
        )
    elif number_of_trafos != 0:
        fig.data = (
            fig.data[start_trafo:start_trafo + number_of_trafos] +
            fig.data[start_trafo + number_of_trafos + 1:]
        )
    elif number_of_lines != 0:
        fig.data = (
            fig.data[:number_of_lines + color_bar] +
            fig.data[number_of_lines + color_bar + 1:]
        )

    # Add the custom arrows to the figure
    fig.add_traces(arrow_traces)

    # --- Final layout adjustments ---
    if on_map:
        map_center_lon, map_center_lat = compute_box_center(bus_coords)
        map_dict = dict(
            center=dict(
                lon=map_center_lon,
                lat=map_center_lat
            ),
            zoom=zoomlevel,
            style=map_style,
        )
        fig.update_layout(map=map_dict)
        config = {}
    else:
        config = {'scrollZoom': True}
        fig.update_layout(width=None, height=None)

    # --- Save the figure ---
    if filename is not None:
        filename = Path(filename)
        if filename.suffix == '.html':
            fig.write_html(filename, config=config)
        elif filename.suffix in ['.svg', '.png', '.jpeg', '.webp', '.pdf']:
            try:
                fig.write_image(filename)  # Requires Kaleido
            except Exception as e:
                warnings.warn(
                    f"Failed to save image: {e}. "
                    "Ensure Kaleido is installed for non-HTML formats."
                )
        else :
            warnings.warn(
                f"Unknown file format for saving Plotly figure: '{filename}'. "
                "Supported formats: '.html', '.svg', '.png', '.jpeg', '.webp', "
                f"'.pdf'."
            )

    return fig