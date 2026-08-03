import pyramses 
import pandas as pd
from pathlib import Path

def create_missing_dirs(file_path: Path):
    """
    Creates all missing directories in the given file path.

    Args:
        file_path (Path): The file path for which to create missing directories.
    """

    if file_path.suffix=='':
        parent_dir = file_path
    else:
        parent_dir = file_path.parent

    parent_dir.mkdir(parents=True, exist_ok=True)


def ramses_extractor_to_df(data:pyramses.extractor,case:pyramses.cfg,
                           time_col='Time',
                           fileName:Path=None)->pd.DataFrame:
    df = pd.DataFrame(columns=['Time'])
    net= case.ppnet
    
    #Get all available signals name from the 'extractor' object's attributes
    branch_list = data._braname
    bus_list    = data._busname
    dctl_list   = data._dctlname
    dctlobs_list= data._dctlobsname
    excobs_list = data._excobsname
    inj_list    = data._injname
    injobs_list = data._injobsname
    #load_list   = data._ldname #Seems to be broken there never something in this list
    shunt_list  = data._shuname
    sync_list   = data._syncname
    torobs_list = data._torobsname
    twop_list   = data._twopname
    twopobs_list= data._twopobsname
    
    #Set the time columns of the dataframe
    t=data._time
    df[time_col]=t
    
    #Renaming convention -- User defined
    conv_bus_name={'mag':'V',
                   'pha':'θ'}
    
    conv_sync_name={'A':'δ',
                    'S':'ω',
                    'FW':'ψf',
                    'DD':'ψd1',
                    'QD':'ψq1',
                    'QW':'ψq2',
                    'FC':'if_mach',
                    'FV':'vf_mach',
                    'T':'Tm_mach',
                    'ET':'Te',
                    'SC':'ω_coi'}
    conv_exc_name={}
    conv_tor_name={}
    
    conv_branch_name={'PF':'P',
                      'PT':'P',
                      'QF':'Q',
                      'QT':'Q',
                      'RM':'n_tfo',
                      'RA':'φ_tfo',
                      }
    conv_inj_name={}
    #conv_load_name={}


    #bus_units={'mag':'pu','pha':'deg'}
    
    for b_n in bus_list:
        b_obj=data.getBus(b_n)
        obs = b_obj.obsdict
        for obs_n in obs:
            attr = getattr(b_obj,obs_n,None)
            if attr is not None:
                obs_split=obs[obs_n].split('(')
                
                if len(obs_split)<=1:
                    unit=''
                else:
                    unit=(obs_split[1]).split(')')[0]
                    
                if obs_n in conv_bus_name.keys():
                    sig_n=conv_bus_name[obs_n]
                else:
                    sig_n=obs_n
                col_name=f'{sig_n}{b_n}_{unit}'
                df[col_name]=attr.value
    
    
    
    for s_n in sync_list:
        s_obj=data.getSync(s_n)
        obs = s_obj.obsdict
        for obs_n in obs:
            attr = getattr(s_obj,obs_n,None)
            if attr is not None:
                obs_split=obs[obs_n].split('(')
                
                if len(obs_split)<=1:
                    unit=''
                else:
                    unit=(obs_split[1]).split(')')[0]
                    
                if obs_n in conv_sync_name.keys():
                    sig_n=conv_sync_name[obs_n]
                else:
                    sig_n=obs_n
                col_name=f'{s_n}.{sig_n}_{unit}'
                df[col_name]=attr.value
                
        exc_obj=data.getExc(s_n)
        obs = exc_obj.obsdict
        for obs_n in obs:
            attr = getattr(exc_obj,obs_n,None)
            if attr is not None:
                obs_split=obs[obs_n].split('(')
                
                if len(obs_split)<=1:
                    unit=''
                else:
                    unit=(obs_split[1]).split(')')[0]
                
                if obs_n in conv_exc_name.keys():
                    sig_n=conv_exc_name[obs_n]
                else:
                    sig_n=obs_n
                col_name=f'{s_n}.{sig_n}_{unit}'
                df[col_name]=attr.value
        
        tor_obj=data.getTor(s_n)
        obs = tor_obj.obsdict
        for obs_n in obs:
            attr = getattr(tor_obj,obs_n,None)
            if attr is not None:
                obs_split=obs[obs_n].split('(')
                
                if len(obs_split)<=1:
                    unit=''
                else:
                    unit=(obs_split[1]).split(')')[0]
                    
                if obs_n in conv_tor_name.keys():
                    sig_n=conv_tor_name[obs_n]
                else:
                    sig_n=obs_n
                col_name=f'{s_n}.{sig_n}_{unit}'
                df[col_name]=attr.value
    
    # for l_n in load_list:
    #     l_obj=data.getLoad(l_n)
    #     obs = l_obj.obsdict
    #     for obs_n in obs:
    #         attr = getattr(l_obj,obs_n,None)
    #         if attr is not None:
    #             unit=(obs[obs_n].split('(')[1]).split(')')[0]
    #             if obs_n in conv_load_name.keys():
    #                 sig_n=conv_load_name[obs_n]
    #             else:
    #                 sig_n=obs_n
    #             col_name=f'{sig_n}_{l_n}_{unit}'
    #             df[col_name]=attr.value
    
    for bra_n in branch_list:
        bra_obj=data.getBranch(bra_n)
        obs = bra_obj.obsdict
        
        
        if (net.line['name']==bra_n).sum()==1:
            from_b = int(net.line.loc[net.line['name']==bra_n,'from_bus'])
            to_b   = int(net.line.loc[net.line['name']==bra_n,'to_bus'])
        elif (net.trafo['name']==bra_n).sum()==1:
            from_b = int(net.trafo.loc[net.trafo['name']==bra_n,'lv_bus'])
            to_b   = int(net.trafo.loc[net.trafo['name']==bra_n,'hv_bus'])
        else:
            raise ValueError(f"Unknow branch name or several branches with"
                             f"the same name, mismatch between pandapower"
                             f" and RAMSES.")
        
        for obs_n in obs:
            attr = getattr(bra_obj,obs_n,None)
            if attr is not None:
                obs_split=obs[obs_n].split('(')
                
                if len(obs_split)<=1:
                    unit=''
                else:
                    unit=(obs_split[1]).split(')')[0]
                
                if obs_n in conv_branch_name.keys():
                    sig_n=conv_branch_name[obs_n]
                else:
                    sig_n=obs_n
                    
                if obs_n.endswith('F'):
                    col_name=f'{sig_n}{from_b}{to_b}_{unit}'
                elif obs_n.endswith('T'):
                    col_name=f'{sig_n}{to_b}{from_b}_{unit}'
                else:
                    col_name=f'{sig_n}{bra_n}_{unit}'
                
                if col_name in df.columns: #mean that this is a second branch in parallel
                    col_name=f'{col_name}_b'
                
                df[col_name]=attr.value
    
    for inj_n in inj_list:
        inj_obj=data.getInj(inj_n)
        obs = inj_obj.obsdict
        
        if inj_n.startswith('L'): #HARCODED clever solution needed
            PQ_factor=100 #MW/MVA/Mvar            
        else:
            PQ_factor=1
        
        for obs_n in obs:
            attr = getattr(inj_obj,obs_n,None)
            if attr is not None:
                obs_split=obs[obs_n].split('(')

                if len(obs_split)<=1:
                    unit=''
                else:
                    unit=(obs_split[1]).split(')')[0]
                
                if obs_n in conv_inj_name.keys():
                    sig_n=conv_inj_name[obs_n]
                else:
                    sig_n=obs_n
                    
                
                col_name=f'{inj_n}.{sig_n}_{unit}'
                
                if obs_n.startswith('P'):
                    col_name=f'{inj_n}.{sig_n}_MW'
                    df[col_name]=attr.value * PQ_factor
                elif obs_n.startswith('Q'):
                    col_name=f'{inj_n}.{sig_n}_Mvar'
                    df[col_name]=attr.value * PQ_factor
                else:
                    df[col_name]=attr.value
                    
    if fileName is not None:
        create_missing_dirs(fileName)    
        df.to_csv(fileName,index=False,float_format='%.6f')

    
    return df