'''
Set all Storage Ring BPMs' settings appropriate for PSD monitoring
'''

import datetime
from cothread.catools import caput, DBR_CHAR_STR

message = "Making BPMs' settings for PSD ..."
print("%s: %s"%(datetime.datetime.now(), str(message)))
#put array of characters as string
caput('SR-APHLA{BPM}PSD:Status-Wf',str(message),datatype=DBR_CHAR_STR)

def _caput(*args):
    caput(*args, repeat_value = True)

#prefix: ID BPMs: SR-HLA{}IDBPMs*; #Arc BPMs: SR-HLA{}AllBPMs*
#9 groups of PVs
machine_type_PVs = ['SR-HLA{}IDBPMs-Loc:Machine-SP', 'SR-HLA{}AllBPMs-Loc:Machine-SP']
trigger_PVs = ['SR-HLA{}IDBPMs-Trig:TrigSrc-SP', 'SR-HLA{}AllBPMs-Trig:TrigSrc-SP']
mode_PVs = ['SR-HLA{}IDBPMs-DDR:WfmSel-SP', 'SR-HLA{}AllBPMs-DDR:WfmSel-SP']
pilot_tone_PVs = ['SR-HLA{}IDBPMs-Rf:PtPwr-SP', 'SR-HLA{}AllBPMs-Rf:PtPwr-SP']
event_code_PVs = ['SR-HLA{}IDBPMs-TimeMode-SP', 'SR-HLA{}AllBPMs-TimeMode-SP']
burst_len_PVs=['SR-HLA{}IDBPMs-Burst:FaEnableLen-SP','SR-HLA{}AllBPMs-Burst:FaEnableLen-SP']
erec_len_PVs=['SR-HLA{}IDBPMs-ERec:FaEnableLen-SP','SR-HLA{}AllBPMs-ERec:FaEnableLen-SP']
ddr_offset_PVs = ['SR-HLA{}IDBPMs-ddrFaOffset', 'SR-HLA{}AllBPMs-ddrFaOffset'] 
switch_mode_PVs = ['SR-HLA{}IDBPMs:SwitchMode-Cmd', 'SR-HLA{}AllBPMs:SwitchMode-Cmd']

pv_lists = [machine_type_PVs, trigger_PVs, mode_PVs, pilot_tone_PVs, event_code_PVs,
            burst_len_PVs, erec_len_PVs, ddr_offset_PVs, switch_mode_PVs]
values = [5, 1, 2, 0, 1, 100100, 100000, 0, 1]
for (pv_list, value) in zip(pv_lists, values):
    _caput(pv_list, value)

caput('SR:C21-PS{Pinger}Mode:Trig-Sel',0)

caput('SR-APHLA{BPM}PSD:IgnoreMatlabTrigger-Cmd',1)  #ignore matlab tirgger
caput('SR-APHLA{BPM}PSD:LiveData-Cmd',1)  #Live data
caput('SR-APHLA{BPM}PSD-Cmd.SCAN',6)  #60 sec scan 



'''
#ID BPMs
MachineTypePV='SR-HLA{}IDBPMs-Loc:Machine-SP'
TriggerPV='SR-HLA{}IDBPMs-Trig:TrigSrc-SP'
modePV='SR-HLA{}IDBPMs-DDR:WfmSel-SP'
PTOnOffPV='SR-HLA{}IDBPMs-Rf:PtPwr-SP'
TBTOffsetPV='SR-HLA{}IDBPMs-ddrTbtOffset'
FAOffsetPV='SR-HLA{}IDBPMs-ddrFaOffset'
TBT_recordPV='SR-HLA{}IDBPMs-ERec:TbtEnableLen-SP'
FA_recordPV='SR-HLA{}IDBPMs-ERec:FaEnableLen-SP'
TBT_burstPV='SR-HLA{}IDBPMs-Burst:TbtEnableLen-SP'
FA_burstPV='SR-HLA{}IDBPMs-Burst:FaEnableLen-SP'
Event_CodePV='SR-HLA{}IDBPMs-TimeMode-SP'

caput(MachineTypePV,5)
caput(TriggerPV,1)
caput(PTOnOffPV,0)
caput(modePV,2)
caput(Event_CodePV,1)


PVlength=['SR-HLA{}IDBPMs-Burst:FaEnableLen-SP','SR-HLA{}IDBPMs-ERec:FaEnableLen-SP','SR-HLA{}IDBPMs-ddrFaOffset']
caput(PVlength,[100100,100000,0])
caput('SR-HLA{}IDBPMs:SwitchMode-Cmd',1)

###Arc BPM

MachineTypePV='SR-HLA{}AllBPMs-Loc:Machine-SP'
TriggerPV='SR-HLA{}AllBPMs-Trig:TrigSrc-SP'
modePV='SR-HLA{}AllBPMs-DDR:WfmSel-SP'
PTOnOffPV='SR-HLA{}AllBPMs-Rf:PtPwr-SP'
TBTOffsetPV='SR-HLA{}AllBPMs-ddrTbtOffset'
FAOffsetPV='SR-HLA{}AllBPMs-ddrFaOffset'
TBT_recordPV='SR-HLA{}AllBPMs-ERec:TbtEnableLen-SP'
FA_recordPV='SR-HLA{}AllBPMs-ERec:FaEnableLen-SP'
TBT_burstPV='SR-HLA{}AllBPMs-Burst:TbtEnableLen-SP'
FA_burstPV='SR-HLA{}AllBPMs-Burst:FaEnableLen-SP'
Event_CodePV='SR-HLA{}AllBPMs-TimeMode-SP'

caput(MachineTypePV,5)
caput(TriggerPV,1)
caput(PTOnOffPV,0)
caput(modePV,2)
caput(Event_CodePV,1)


PVlength=['SR-HLA{}AllBPMs-Burst:FaEnableLen-SP','SR-HLA{}AllBPMs-ERec:FaEnableLen-SP','SR-HLA{}AllBPMs-ddrFaOffset']
caput(PVlength,[100100,100000,0])
caput('SR-HLA{}AllBPMs:SwitchMode-Cmd',1)
'''
