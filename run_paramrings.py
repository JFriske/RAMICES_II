import os
import subprocess
import numpy as np
from datetime import datetime


beg = datetime.now()
#readout time at beginning and output at the end



fhnsmarr = np.linspace(0.25, 0.9, 10)
fhccsnarr = 1.0 -(0.5*(1.0 - fhnsmarr))
fhagmarr = np.linspace(0.1,0.65, 10)
fhsn1aarr = (0.85, 0.9, 0.93, 0.95, 0.97, 0.98, 0.99, 0.995, 0.999, 0.9999)
ejectglobalarr = np.linspace(0.3,0.6,10)
ejectnucleararr = np.linspace(0.45,0.95,10)

inflowarr = [0,1]

paramarr = [fhccsnarr, fhnsmarr, fhagmarr, fhsn1aarr, ejectglobalarr, ejectnucleararr]
defaultarr = [0.8, 0.6, 0.7, 0.99, 0.55, 0.65]

total_param = 6
arrlength = 10

fullparamarr = np.empty(0)

for paramnr, param in enumerate(paramarr):
    dummy_default = np.full(arrlength, defaultarr[paramnr])
    # if paramnr == 0:
    #     sequence = np.append(paramarr[paramnr],np.tile(dummy_default,arrlength-paramnr-1) )
    # elif paramnr ==9:
    #     sequence = np.append(np.tile(dummy_default,paramnr),paramarr[paramnr])
    # else:    
    sequence = np.append(np.tile(dummy_default,paramnr),paramarr[paramnr])
    sequence =  np.append(sequence, np.tile(dummy_default,total_param-paramnr-1) )
    print(sequence)
    fullparamarr = np.append(fullparamarr, sequence, axis = 0 )
    #print(fullparamarr[paramnr]

fullparamarr = np.reshape(fullparamarr, (-1,total_param*arrlength))    
print(fullparamarr)
#print(fullparamarr.shape())
print(len(fullparamarr))

suitename = "NucParamRingTest"

outputdir = "/disk/xray8/jksf/ChemicalEvolution/"
#launchmake = subprocess.Popen("../../make", shell=True, stdout=subprocess.PIPE)
#launchmake.wait()

nrglob = 0
nrnuc = 0 
processes = []
logfilearr = []
for fhccsn, fhnsm, fhagb, fhsn1a, ejectglob, ejectnuc in zip(fullparamarr[0], fullparamarr[1], fullparamarr[2], fullparamarr[3], fullparamarr[4], fullparamarr[5]):
    filenameglob = f"Global_fhccsn{fhccsn:.3f}_fhnsm{fhnsm:.3f}_fhagb{fhagb:.3f}_fhsn1a{fhsn1a:.4f}_ejectglob{ejectglob:.3f}_ejectnuc{ejectnuc:.3f}"
    with open ("config/"+ suitename +"/" + filenameglob+".config", "w") as outfile:
        with open ("config/"+ suitename +"/stdglob_paramtest.config") as infile:
            outfile.write("output "+outputdir+ "Output/"+suitename + "/" +filenameglob+"\n")
            outfile.write(infile.read())
            outfile.write("\n")
            outfile.write("fh-ccsn {:.4f}\n".format(fhccsn))
            outfile.write("fh-nsm {:.4f}\n".format(fhnsm))
            outfile.write("fh-agb {:.4f}\n".format(fhagb))
            outfile.write("fh-sn1a {:.4f}\n".format(fhsn1a))
            outfile.write("eject {:.4f}\n".format(ejectglob))
                            

    logfile = outputdir+ "Output/"+suitename + "/" +filenameglob+ "/output.log"

    outputfolder = outputdir + "/" + "Output/"+suitename+ "/" + filenameglob
    if not os.path.exists(outputfolder):
        os.mkdir(outputfolder)
    launchglob = subprocess.Popen("./Ramices_Launch.sh -config config/"+suitename + "/" + filenameglob + ".config", shell=True, stdout= subprocess.PIPE, bufsize=1)
    processes.append(launchglob)
    logfilearr.append(logfile)



for  nrglob, (launchglob, logfile) in enumerate(zip(processes, logfilearr)):
    with open(logfile, 'w') as log:
        for line in iter(launchglob.stdout.readline, b''):
            #print ("glob" + str(nrglob) + " " + str(line))
            log.write(str(line))
    launchglob.stdout.close()
    launchglob.wait()
    print("Done Global Nr " + str(nrglob) + "\n")

print(len(processes))
processes = []
logfilearr = []
print(len(processes))

for fhccsn, fhnsm, fhagb, fhsn1a, ejectglob, ejectnuc in zip(fullparamarr[0], fullparamarr[1], fullparamarr[2], fullparamarr[3], fullparamarr[4], fullparamarr[5]):
    for inflow in inflowarr:
        filenamenuc = f"Nuclear_fhccsn{fhccsn:.3f}_fhnsm{fhnsm:.3f}_fhagb{fhagb:.3f}_fhsn1a{fhsn1a:.4f}_ejectglob{ejectglob:.3f}_ejectnuc{ejectnuc:.3f}"
        if inflow == 1:
            filenamenuc = f"Nuclear_fhccsn{fhccsn:.3f}_fhnsm{fhnsm:.3f}_fhagb{fhagb:.3f}_fhsn1a{fhsn1a:.4f}_ejectglob{ejectglob:.3f}_ejectnuc{ejectnuc:.3f}_inflow"
        filenameglob = f"Global_fhccsn{fhccsn:.3f}_fhnsm{fhnsm:.3f}_fhagb{fhagb:.3f}_fhsn1a{fhsn1a:.4f}_ejectglob{ejectglob:.3f}_ejectnuc{ejectnuc:.3f}"
        with open ("config/"+ suitename +"/" + filenamenuc+".config", "w") as outfile:
            with open ("config/"+ suitename +"/stdnuc_paramtest.config") as infile:
                outfile.write("output "+outputdir+ "/Output/"+suitename + "/" +filenamenuc+"\n")
                outfile.write("readin-dir "+outputdir+ "/Output/"+suitename + "/"  +filenameglob + "\n")
                outfile.write(infile.read())
                outfile.write("\n")
                outfile.write("fh-ccsn {:.4f}\n".format(fhccsn))
                outfile.write("fh-nsm {:.4f}\n".format(fhnsm))
                outfile.write("fh-agb {:.4f}\n".format(fhagb))
                outfile.write("fh-sn1a {:.4f}\n".format(fhsn1a))
                outfile.write("eject {:.4f}\n".format(ejectnuc))
                outfile.write("inflow-on {}\n".format(inflow))

        logfile = outputdir+ "Output/"+suitename + "/" +filenamenuc+ "/output.log"

        outputfolder = outputdir + "/" + "Output/"+suitename+ "/" + filenamenuc
        if not os.path.exists(outputfolder):
            os.mkdir(outputfolder)
        launchnuc = subprocess.Popen("./Ramices_Launch.sh -config config/"+suitename + "/" + filenamenuc + ".config", shell=True, stdout= subprocess.PIPE, bufsize=1)
        processes.append(launchnuc)
        logfilearr.append(logfile)



for  nrnuc, (launchnuc, logfile) in enumerate(zip(processes, logfilearr)):
    with open(logfile, 'w') as log:
        for line in iter(launchnuc.stdout.readline, b''):
            #print ("glob" + str(nrglob) + " " + str(line))
            log.write(str(line))
    launchnuc.stdout.close()
    launchnuc.wait()
    print("Done Nuclear Nr " + str(nrnuc) + "\n")

                                
end = datetime.now() -beg

current_time = end#.strftime("%H:%M:%S")
print("Current Time =", current_time)