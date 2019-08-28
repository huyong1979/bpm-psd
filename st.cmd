#!/epics/iocs/elauncher/elauncher/bin/linux-x86_64/scriptlaunch

epicsEnvSet("EPICS_CA_ADDR_LIST", "10.0.153.255")
epicsEnvSet("EPICS_CA_AUTO_ADDR_LIST", "NO")
epicsEnvSet("PATH", "/home/skongtawong/anaconda2/bin:$PATH")
#epicsEnvSet("PYTHONPATH", "/epics/op/apps/aphla/lib/python")
#epicsEnvSet("APHLA_CONFIG_DIR", "/epics/data/aphla/apconf")

#cd "/epics/iocs/elauncher/elauncher"
cd "/epics/iocs/bpm-psd"
dbLoadDatabase("./scriptlaunch.dbd",0,0)
scriptlaunch_registerRecordDeviceDriver(pdbbase)

dbLoadRecords("/usr/lib/epics/db/iocAdminSoft.db", "IOC=SR-APHLA{IOC:BPMPSD}")
dbLoadRecords ("/usr/lib/epics/db/save_restoreStatus.db", "P=SR-APHLA{IOC:BPMPSD}")
dbLoadTemplate("psd.substitutions")
dbLoadRecords("psdMisc.db")

save_restoreSet_status_prefix("SR-APHLA{IOC:BPMPSD}")
set_savefile_path("./as", "/save")
set_requestfile_path("./as", "/req")
set_pass0_restoreFile("settings_pass0.sav")
set_pass1_restoreFile("settings_pass1.sav")

#asSetFilename("/cf-update/acf/default.acf")

iocInit()

makeAutosaveFileFromDbInfo("./as/req/settings_pass0.req", "autosaveFields_pass0")
create_monitor_set("settings_pass0.req", 30 , "")
makeAutosaveFileFromDbInfo("./as/req/settings_pass1.req", "autosaveFields_pass1")
create_monitor_set("settings_pass1.req", 30 , "")

#caPutLogInit("ioclog.cs.nsls2.local:7004", 1)
