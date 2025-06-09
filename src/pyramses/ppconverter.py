import pandapower as pp
import os 
from pathlib import Path
import warnings

from .globals import RAMSESError, __runTimeObs__, CustomWarning, silentremove

warnings.showwarning = CustomWarning

def convertBusDatatoPP(net, bus_list_of_str:list):
    print(f"Bus list: {bus_list_of_str}")
    print(f"Number of buses: {len(bus_list_of_str)}")
    print(f"Type of bus_list_of_str: {type(bus_list_of_str)}")
    for bus_str in bus_list_of_str:
        # Split the bus line into components, the separator is a space or mutiple spaces
        bus_components = bus_str.split()
        print(f"Bus components: {bus_components}")

        if len(bus_components) == 2:
            bus_name = bus_components[0].strip()
            vn_kv = float(bus_components[1].strip())
            pload_mw = 0.0
            qload_mvar = 0.0
            bshunt_mvar = 0.0
            qshunt_mvar = 0.0
        elif len(bus_components) == 6:
            bus_name = bus_components[0].strip()
            vn_kv = float(bus_components[1].strip())
            pload_mw    = float(bus_components[2].strip())
            qload_mvar  = float(bus_components[3].strip())
            bshunt_mvar = float(bus_components[4].strip())
            qshunt_mvar = float(bus_components[5].strip())
        else:
            raise ValueError(f"Invalid bus definition: {bus_str}. Expected 2 or 6 components.")
        
        qload_tot_mvar = qload_mvar + qshunt_mvar 
        #B_shunt is discarded in pandapower, there is something similar with
        # "const_z_percent" and "const_i_percent" in pandapower, but it is not
        #  not used in the same way as in PFC of Stepss.

        bus_id = pp.create_bus(net, vn_kv, name=bus_name,index=int(bus_name))
        pp.create_load(net, bus_id, pload_mw, q_mvar=qload_tot_mvar,
                        name='load_'+bus_name)
    #debug : print network buses and loads
    print(f"Created buses in pandapower network:\n {net.bus}")
    print(f"Created loads in pandapower network:\n {net.load}")

def convertLineDatatoPP(net, line_list_of_str:list, length_km_list=None):
    print(f"Line list: {line_list_of_str}")
    print(f"Number of lines: {len(line_list_of_str)}")
    print(f"Type of line_list_of_str: {type(line_list_of_str)}")
    
    #get network nominal frequency :
    f_nom = net.f_hz if hasattr(net, 'f_hz') else 50.0  # Default to 50.0 Hz if not set

    l=0
    for line_str in line_list_of_str:
        # Split the line line into components, the separator is a space or mutiple spaces
        line_components = line_str.split()
        print(f"Line components: {line_components}")

        if len(line_components) == 8:
            line_name = line_components[0].strip()
            from_bus_name = line_components[1].strip()
            to_bus_name = line_components[2].strip()
            r_ohm = float(line_components[3].strip())
            x_ohm = float(line_components[4].strip())
            wc2_uS = float(line_components[5].strip())
            sn_l_mva = float(line_components[6].strip())
            br_status = bool(int(line_components[7].strip()))  # Convert to True/False
            
        else:
            raise ValueError(f"Invalid line definition: {line_str}. Expected 8 components.")
        
        #Create the line in pandapower
        from_bus_id = net.bus[net.bus.name == from_bus_name].index[0]
        to_bus_id = net.bus[net.bus.name == to_bus_name].index[0]
        
        #Maximum current in kA
        bus_voltage_kv = net.bus.vn_kv[from_bus_id] # Assuming both buses have the same voltage level
        max_i_ka = sn_l_mva / (bus_voltage_kv)      # Convert MVA to kA using the formula I = S / V

        if length_km_list is not None:
            length_km = length_km_list[l] if l < len(length_km_list) else 1.0  # Default to 1.0 km if not provided
            l+=1
        else:
            length_km = 1.0 # Default to 1.0 km if not provided

    
        r_ohm_per_km = r_ohm / length_km
        x_ohm_per_km = x_ohm / length_km
        c_nf_per_km = (( (2* wc2_uS * 1e-6)/ (f_nom * 2 * 3.14159) )/ length_km) /1e-9 # Convert wc2_uS to c_nf_per_km  
        pp.create_line_from_parameters(net, from_bus_id, to_bus_id, length_km,
                                        r_ohm_per_km, x_ohm_per_km, c_nf_per_km,
                                            max_i_ka, name=line_name,  in_service=br_status)
    #debug : print network lines
    print(f"Created lines in pandapower network:\n {net.line}")

def createTfo(net,from_bus_name:str, to_bus_name:str, r_from_per_base:float,
               x_from_per_base:float, b_to_per_base:float, n_per:float,
                 sn_mva:float, shift_degree:float=0.0, tap_side=None,tap_neutral=None,
                 tap_max=None, tap_min=None, tap_step_percent=None, transfo_name=None,in_service=True):
    #get network nominal frequency :
    f_nom = net.f_hz if hasattr(net, 'f_hz') else 50.0  # Default to 50.0 Hz if not set
    sn_mva_net = net.sn_mva if hasattr(net, 'sn_mva') else 100.0  # Default to 100.0 MVA if not set

    # According to Stepss doc and Pandapower,
    # FROM bus = Low Voltage (LV) bus, STEPSS/Pandapower
    # TO bus = High Voltage (HV) bus, STEPSS/Pandapower

    #Create the transformer in pandapower
    lv_bus_id = net.bus[net.bus.name == from_bus_name].index[0]
    hv_bus_id = net.bus[net.bus.name == to_bus_name].index[0]

    vb_kv_hv = net.bus.vn_kv[hv_bus_id]  # High voltage bus voltage
    vb_kv_lv = net.bus.vn_kv[lv_bus_id]  # Low voltage bus voltage

    voc_kv_lv= vb_kv_lv
    voc_kv_hv= n_per/100 * (voc_kv_lv/vb_kv_lv) * vb_kv_hv  # Voltage at the HV side in kV
    print(f"Voltage at HV side: {voc_kv_hv} kV")
    print(f"Voltage at LV side: {voc_kv_lv} kV")
    print(f"Ratio of transformer: {n_per}, from {vb_kv_lv} kV to {vb_kv_hv} kV")


    r_pu = r_from_per_base /100 *(vb_kv_lv/voc_kv_lv)**2
    x_pu = x_from_per_base /100 *(vb_kv_lv/voc_kv_lv)**2 #From Stepss documentation, R and X are the resistances and reactances at the low voltage side, in pu
    
    p_mw_cu = r_pu *sn_mva
    vkr_percent = (p_mw_cu / sn_mva) * 100  # Convert to percentage

    z_pu = complex(r_pu, x_pu)  # Impedance in pu
    z_mag_pu = abs(z_pu)  # Magnitude of the impedance in pu

    vk_percent = z_mag_pu * 100  * (sn_mva / sn_mva_net)  # Convert to percentage

    b_pu = b_to_per_base / 100 * (vb_kv_lv/voc_kv_lv)**2  #From Stepss documentation, B is the shunt admittance at the low voltage side
    y_pu = complex(0, b_pu)  # Admittance in pu
    y_mag_pu = abs(y_pu)
    i0_percent = 100 * y_mag_pu

    pfe_kw = 0.0  # Power loss in kW, Conductance neglected in STEPSS, set to 0

    

    pp.create_transformer_from_parameters(net, hv_bus_id, lv_bus_id,
                                            sn_mva, voc_kv_hv, voc_kv_lv,
                                            vkr_percent, vk_percent, pfe_kw,
                                            i0_percent, shift_degree=shift_degree,
                                            tap_side=tap_side, tap_neutral=tap_neutral,
                                            tap_max=tap_max, tap_min=tap_min,
                                            tap_step_percent=tap_step_percent,
                                            name=transfo_name, in_service=in_service,
                                            )


def convertTransfoDatatoPP(net, transfo_list_of_str:list, trfo_list_of_str:list, pshift_list_of_str:list=None, ltc_v_list_of_str:list=None):
    
    for transfo_str in transfo_list_of_str:
        # Split the transformer line into components, the separator is a space or mutiple spaces
        transfo_components = transfo_str.split()
        print(f"Transformer components: {transfo_components}")

        if len(transfo_components) == 11:
            #TRANSFO NAME FROMBUS TOBUS R X B1 B2 N SNOM SHIFT BR ;
            #Pi model transformer
            trafo_name = transfo_components[0].strip()
            from_bus_name = transfo_components[1].strip() 
            to_bus_name = transfo_components[2].strip()
            r_from_per_base = float(transfo_components[3].strip())
            x_from_per_base = float(transfo_components[4].strip())
            b_from_per_base = float(transfo_components[5].strip())#Discarded in pandapower
            b_to_per_base = float(transfo_components[6].strip())
            n_per = float(transfo_components[7].strip())
            sn_mva = float(transfo_components[8].strip())
            shift_degree = float(transfo_components[9].strip())
            br_status = bool(int(transfo_components[10].strip()))  # Convert to True/False
           
        else:
            raise ValueError(f"Invalid transformer definition: {transfo_str}. Expected 11 components.")
        
        createTfo(net, from_bus_name, to_bus_name, r_from_per_base,
                   x_from_per_base, b_to_per_base, n_per, sn_mva,
                   shift_degree=shift_degree, transfo_name=trafo_name,in_service=br_status)
        
    for trfo_str in trfo_list_of_str:
        # Split the transformer line into components, the separator is a space or mutiple spaces
        trfo_components = trfo_str.split()
        print(f"Transformer components: {trfo_components}")

        if len(trfo_components) == 15:
            #TRFO NAME FROMBUS TOBUS CONBUS R X B N SNOM NFIRST NLAST NBPOS TOLV VDES BR ;
            #Pi model transformer with tap changer
            transfo_name = trfo_components[0].strip()
            from_bus_name = trfo_components[1].strip()
            to_bus_name = trfo_components[2].strip()
            con_bus_name = trfo_components[3].strip() #not used in pandapower, used by ramses controller
            r_from_per_base = float(trfo_components[4].strip())
            x_from_per_base = float(trfo_components[5].strip())
            b_to_per_base = float(trfo_components[6].strip())
            n_per = float(trfo_components[7].strip())
            sn_mva = float(trfo_components[8].strip())
            n_first_per = float(trfo_components[9].strip())
            n_last_per = float(trfo_components[10].strip())
            nb_pos = int(trfo_components[11].strip())
            tol_v_pu = float(trfo_components[12].strip()) #not used in pandapower, used by ramses controller
            v_des_pu = float(trfo_components[13].strip()) #not used in pandapower, used by ramses controller
            br_status = bool(int(trfo_components[14].strip()))  # Convert to True/False
            shift_degree =0.0  # No phase shift for tap changers in pandapower
            
            #Compute Tap position of the trfo
            p_neutral=int((n_per - n_first_per)*(nb_pos-1)/(n_last_per-n_first_per)+1) 
            print (f"Tap position: {p_neutral} for transformer {trafo_name}")    
            tap_side = 'lv'

            tap_min = 1
            tap_step_percent = (n_last_per - n_first_per) / (nb_pos-1)  # Tap step in percentage
            tap_max = int(tap_min + ((nb_pos-1) * tap_step_percent))  # Maximum tap position
        else: 
            raise ValueError(f"Invalid transformer definition: {trfo_str}. Expected 15 components.")
        
        createTfo(net, from_bus_name, to_bus_name, r_from_per_base,
                   x_from_per_base, b_to_per_base, n_per, sn_mva,
                   shift_degree=shift_degree, tap_side=tap_side, tap_neutral=p_neutral,
                   tap_max=tap_max, tap_min=tap_min, tap_step_percent=tap_step_percent,
                   transfo_name=transfo_name, in_service=br_status)
        
    for pshift_str in pshift_list_of_str:
        print(f'Not implemented yet:\n PSHIFT-P {pshift_str}')
        break

    for ltc_v_str in ltc_v_list_of_str:
        #LTC-V NAME CONBUS NFIRST NLAST NBPOS TOLV VDES ;
        ltc_v_components = ltc_v_str.split()
        print(f"Load Tap Changer components: {ltc_v_components}")
        if len(ltc_v_components) == 7 or len(ltc_v_components) == 9:
            ltc_name = ltc_v_components[0].strip()
            con_bus_name = ltc_v_components[1].strip() #not used in pandapower, used by ramses controller
            n_first_per = float(ltc_v_components[2].strip())
            n_last_per = float(ltc_v_components[3].strip())
            nb_pos = int(ltc_v_components[4].strip())
            tol_v_pu = float(ltc_v_components[5].strip())# not used in pandapower, used by ramses controller
            v_des_pu = float(ltc_v_components[6].strip())# not used in pandapower, used by ramses controller

            #Get the pp transformer characteristics 
            transfo_row = net.trafo[net.trafo.name == ltc_name]
            print(f"Transformer row: {transfo_row}")
            print(f"type of transfo_row: {type(transfo_row)}")
            if transfo_row.empty:
                raise ValueError(f"Transformer {ltc_name} not found in pandapower network.")
            voc_kv_hv = transfo_row['vn_hv_kv'].values[0]  # High voltage side voltage
            voc_kv_lv = transfo_row['vn_lv_kv'].values[0]  # Low voltage side voltage
            n_per = (voc_kv_hv / voc_kv_lv) * 100  # Transformer ratio in percentage 
            print(n_per)
            print(type(n_per))
            #Compute Tap position of the ltc
            p_neutral=int((n_per - n_first_per)*(nb_pos-1)/(n_last_per-n_first_per)+1) 
            print (f"Tap position: {p_neutral} for transformer {trafo_name}")    
            tap_side = 'lv'

            tap_min = 1
            tap_step_percent = (n_last_per - n_first_per) / (nb_pos-1)  # Tap step in percentage
            tap_max = int(tap_min + ((nb_pos-1) * tap_step_percent))  # Maximum tap position

            #Update the transformer parameters in pandapower
            net.trafo.loc[net.trafo.name == ltc_name, 'tap_side'] = tap_side
            net.trafo.loc[net.trafo.name == ltc_name, 'tap_neutral'] = p_neutral
            net.trafo.loc[net.trafo.name == ltc_name, 'tap_max'] = tap_max
            net.trafo.loc[net.trafo.name == ltc_name, 'tap_min'] = tap_min
            net.trafo.loc[net.trafo.name == ltc_name, 'tap_step_percent'] = tap_step_percent

        else:
            raise ValueError(f"Invalid load tap changer definition: {ltc_v_str}. Expected 7 components.")
        
    #debug : print network transformers
    print(f"Created transformers in pandapower network:\n {net.trafo}")   

def convertDataToPandaPowerNetwork(datfiles_list:list,net_name='pyramses_network'):
    """
    Converts a list of .dat files to a pandapower network.

    Parameters
    ----------
    datfiles_list : list of str
        List of paths to .dat files.

    Returns
    -------
    pandapowerNet : pandapowerNet
        A pandapower network object.
    """

    if not isinstance(datfiles_list, list):
        raise TypeError("datfiles_list must be a list of strings representing file paths")
    
    #Initialize an empty string to hold the data
    data_str = ""

    for datfile in datfiles_list:
        if not datfile.endswith('.dat'):
            raise ValueError(f"File {datfile} is not a .dat file")
        
        if not Path(datfile).resolve().is_file():
            raise FileNotFoundError(f"File {datfile} does not exist")
        
        #Read the file
        with open(datfile, 'r') as file:
            data_str += file.read() #string with all the text from the file
    print (f"Data from {datfile}:\n{data_str}\n")
    #print(f"Type of reading{type(data_str)}\n")
    
    #Remove lines of comments that start with # and solvers options that start with $
    data_lines = data_str.splitlines()
    data_lines = [line for line in data_lines if not line.startswith('#')]
    data_lines = [line for line in data_lines if not line.startswith('$')]

    #Join the lines back into a single string
    data_str = "\n".join(data_lines)

    #Split the string into sentences that end with a semicolon but that can
    #span multiple lines.
    data_str = data_str.replace('\n', '')  # remove newlines
    data_str = data_str.replace(';', ';\n')  # Replace semicolon with semicolon and newline

    print(f"Data after removing comments and solver options:\n{data_str}\n")

    #Split the string into lines to get each sentence
    data_lines=data_str.splitlines() 

    #Starting Keywords separation
    keywords = ['FNOM','BUS', 'LINE', 'SWITCH', 'TRANSFO', 'NRTP', 'SHUNT',
                 'GENER', 'SLACK', 'SVC', 'LTC-V','PSHIFT-P'] 
    # keywords from Stepss documentation, section PFC only
    
    data_dict = {keyword: [] for keyword in keywords}

    #Assign each line to the corresponding keyword
    for line in data_lines:
        line = line.strip()  # Remove leading and trailing whitespace
        if not line:  # Skip empty lines
            continue
        for keyword in keywords:
            if line.startswith(keyword):
                #remove the keyword from the line
                line = line.replace(keyword, '').strip()
                #remove the semicolon at the end of the line
                line = line.replace(';', '').strip()

                data_dict[keyword].append(line)
                break  # Stop checking after the first match

    print("Data dictionary with separated keywords:")
    for keyword, lines in data_dict.items():
        print(f"{keyword}: {lines}")            

    #pu base of pp :
    sn_pp_mva = 100.0  # Default base power in MVA for pandapower

    #Get grid nominal frequency :
    f_nom_str = (data_dict.get('FNOM', ["50"])[0]).strip()  # Default to 50 if not found
    f_nom = float(f_nom_str) if f_nom_str else 50.0  # Convert to float, default to 50.0 if empty



    net = pp.create_empty_network(name=net_name, f_hz=f_nom, sn_mva=sn_pp_mva)

    # Create buses
    bus_list = data_dict.get('BUS', [])
    convertBusDatatoPP(net, bus_list)
    
    #Create lines
    line_list = data_dict.get('LINE', []) 
    convertLineDatatoPP(net, line_list)

    #Create transformers
    transfo_list = data_dict.get('TRANSFO', [])
    trfo_list = data_dict.get('TRFO', [])
    pshift_list = data_dict.get('PSHIFT-P', [])
    ltc_v_list = data_dict.get('LTC-V', [])
    convertTransfoDatatoPP(net, transfo_list,trfo_list, pshift_list, ltc_v_list)
    

    # #get availble std types for lines
    # std_types = pp.available_std_types(net, element='line') #return a table
    # print(f"Available standard types for lines:\n {std_types}")

    #Create Transformers       

    return net