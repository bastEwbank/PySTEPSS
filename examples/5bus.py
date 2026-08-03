import pyramses
#import pandapower as pp
from pathlib import Path

###CONFIGURE THE CASE STUDY HERE
case_name  = "5bus"
case_dir   = Path().cwd()/case_name 

data_files = ['lf2',
              'solveroptions',
              'dyn',
             ] 
data_files =[f"{f_name}.dat" for f_name in data_files]

obs_file         = "obs.dat"
init_name        = "init.trace"
dist_name        = "events.dst"
trajectory_name  = "out.trj"
cont_trace_name  = "cont.trace"
disc_trace_name  = "disc.trace"
out_trace_name   = "output.trace"
cmd_txt_name     = "cmd.txt"
stepss_cfg_name  = "sim.cfg"
lf_res_name      = "lf_res.dat"

###INITIALIZE THE CASE OBJECT
case = pyramses.cfg()
for file in data_files:
    case.addData(str(case_dir/file))
    
case.addObs(str(case_dir  / obs_file))
case.addInit(str(case_dir / init_name))
case.addDst(str(case_dir  / dist_name))
case.addTrj(str(case_dir  / trajectory_name))
case.addCont(str(case_dir / cont_trace_name))
case.addDisc(str(case_dir / disc_trace_name))
case.addOut(str(case_dir / out_trace_name))

#(Optionnal) Save the case to a file for later import through Python PyRamses script
case.writeCmdFile(str(case_dir/cmd_txt_name))

#(Optionnal)Save the case configuration file, for later import in STEPSS java GUI
case.writeStepssCfgFile(str(case_dir/stepss_cfg_name))

###CONVERT THE CASE TO PANDAPOWER NETWORK and RUN THE POWER FLOW
case.UpdatePandaPowerNetwork(name=case_name)
net=case.getPandaPowerNetwork()
net=case.runPFC(add_results=True,
                tmpFile_path_str=str(case_dir / lf_res_name),
                max_iteration=1000)

###RUN DYNAMIC SIMULATION
sim=pyramses.sim()
ret=sim.execSim(case)
print(ret)

### EXTRACT RESULTS FROM FILE TO NUMPY ARRAY OR DATAFRAME
data = pyramses.extractor(case.getTrj())
print(data)