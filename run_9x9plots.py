import os
import subprocess
import numpy as np
from datetime import datetime


beg = datetime.now()
#readout time at beginning and output at the end



fhccsnarr = np.linspace(0.25, 0.9, 9)
fhnsmarr = 1.0 -(2.0*(1.0 - fhccsnarr))
fhagbarr = np.linspace(0.1,0.65, 9)
fhsn1aarr = (0.85, 0.9, 0.93, 0.95, 0.97, 0.99, 0.995, 0.999, 0.9999)
ejectglobalarr = np.linspace(0.3,0.6,9)

inflowarr = [0,1]

#paramarr = [fhccsnarr, fhnsmarr, fhagmarr, fhsn1aarr, ejectglobalarr, ejectnucleararr]
defaultarr = [0.8, 0.6, 0.7, 0.99, 0.55, 0.65]

nrglob = 0
nrnuc = 0 
processes = []
logfilearr = []

suitename = "NucParam9x9plots"

outputdir = "/disk/xray8/jksf/ChemicalEvolution/"



for fhccsn in fhccsnarr:
    for fhnsm in fhnsmarr:
        filenameglob = f"Global_fhccsn{fhccsn:.3f}_fhnsm{fhnsm:.3f}"
        with open ("config/"+ suitename +"/" + filenameglob+".config", "w") as outfile:
            with open ("config/"+ suitename +"/stdglob_paramtest.config") as infile:
                outfile.write("output "+outputdir+ "Output/"+suitename + "/" +filenameglob+"\n")
                outfile.write(infile.read())
                outfile.write("\n")
                outfile.write("fh-ccsn {:.4f}\n".format(fhccsn))
                outfile.write("fh-nsm {:.4f}\n".format(fhnsm))
                outfile.write("fh-agb {:.4f}\n".format(defaultarr[2]))
                outfile.write("fh-sn1a {:.4f}\n".format(defaultarr[3]))
                outfile.write("eject {:.4f}\n".format(defaultarr[4]))
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


for fhccsn in fhccsnarr:
    for fhagb in fhagbarr:
        filenameglob = f"Global_fhccsn{fhccsn:.3f}_fhagb{fhagb:.3f}"
        with open ("config/"+ suitename +"/" + filenameglob+".config", "w") as outfile:
            with open ("config/"+ suitename +"/stdglob_paramtest.config") as infile:
                outfile.write("output "+outputdir+ "Output/"+suitename + "/" +filenameglob+"\n")
                outfile.write(infile.read())
                outfile.write("\n")
                outfile.write("fh-ccsn {:.4f}\n".format(fhccsn))
                outfile.write("fh-nsm {:.4f}\n".format(defaultarr[1]))
                outfile.write("fh-agb {:.4f}\n".format(fhagb))
                outfile.write("fh-sn1a {:.4f}\n".format(defaultarr[3]))
                outfile.write("eject {:.4f}\n".format(defaultarr[4]))
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

for fhccsn in fhccsnarr:
    for fhsn1a in fhsn1aarr:
        filenameglob = f"Global_fhccsn{fhccsn:.3f}_fhsn1a{fhsn1a:.4f}"
        with open ("config/"+ suitename +"/" + filenameglob+".config", "w") as outfile:
            with open ("config/"+ suitename +"/stdglob_paramtest.config") as infile:
                outfile.write("output "+outputdir+ "Output/"+suitename + "/" +filenameglob+"\n")
                outfile.write(infile.read())
                outfile.write("\n")
                outfile.write("fh-ccsn {:.4f}\n".format(fhccsn))
                outfile.write("fh-nsm {:.4f}\n".format(defaultarr[1]))
                outfile.write("fh-agb {:.4f}\n".format(defaultarr[2]))
                outfile.write("fh-sn1a {:.4f}\n".format(fhsn1a))
                outfile.write("eject {:.4f}\n".format(defaultarr[4]))
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

for fhccsn in fhccsnarr:
    for ejectglob in ejectglobalarr:
        filenameglob = f"Global_fhccsn{fhccsn:.3f}_ejectglob{ejectglob:.3f}"
        with open ("config/"+ suitename +"/" + filenameglob+".config", "w") as outfile:
            with open ("config/"+ suitename +"/stdglob_paramtest.config") as infile:
                outfile.write("output "+outputdir+ "Output/"+suitename + "/" +filenameglob+"\n")
                outfile.write(infile.read())
                outfile.write("\n")
                outfile.write("fh-ccsn {:.4f}\n".format(fhccsn))
                outfile.write("fh-nsm {:.4f}\n".format(defaultarr[1]))
                outfile.write("fh-agb {:.4f}\n".format(defaultarr[2]))
                outfile.write("fh-sn1a {:.4f}\n".format(defaultarr[3]))
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

for fhccsn in fhccsnarr:
    for fhnsm in fhnsmarr:
        for inflow in inflowarr:
            filenamenuc = f"Nuclear_fhccsn{fhccsn:.3f}_fhnsm{fhnsm:.3f}"
            if inflow == 1:
                filenamenuc = f"Nuclear_fhccsn{fhccsn:.3f}_fhnsm{fhnsm:.3f}_inflow"
            filenameglob = f"Global_fhccsn{fhccsn:.3f}_fhnsm{fhnsm:.3f}"
            outputfolder = outputdir + "/" + "Output/"+suitename+ "/" + filenamenuc
            if not os.path.exists(outputfolder):
                os.mkdir(outputfolder)
            with open ("config/"+ suitename +"/" + filenamenuc+".config", "w") as outfile:
                with open ("config/"+ suitename +"/stdnuc_paramtest.config") as infile:
                    outfile.write("output "+outputdir+ "Output/"+suitename + "/" +filenamenuc+"\n")
                    outfile.write("readin-dir "+outputdir+ "/Output/"+suitename + "/"  +filenameglob + "\n")
                    outfile.write(infile.read())
                    outfile.write("\n")
                    outfile.write("fh-ccsn {:.4f}\n".format(fhccsn))
                    outfile.write("fh-nsm {:.4f}\n".format(fhnsm))
                    outfile.write("fh-agb {:.4f}\n".format(defaultarr[2]))
                    outfile.write("fh-sn1a {:.4f}\n".format(defaultarr[3]))
                    outfile.write("eject {:.4f}\n".format(defaultarr[5]))
                    outfile.write("inflow-on {}\n".format(inflow))
            logfile = outputdir+ "Output/"+suitename + "/" +filenamenuc+ "/output.log"

            
            launchglob = subprocess.Popen("./Ramices_Launch.sh -config config/"+suitename + "/" + filenamenuc + ".config", shell=True, stdout= subprocess.PIPE, bufsize=1)
            processes.append(launchglob)
            logfilearr.append(logfile)


for  nrnuc, (launchnuc, logfile) in enumerate(zip(processes, logfilearr)):
    with open(logfile, 'w') as log:
        for line in iter(launchnuc.stdout.readline, b''):
            #print ("glob" + str(nrglob) + " " + str(line))
            log.write(str(line))
    launchnuc.stdout.close()
    launchnuc.wait()
    print("Done Nuclear Nr " + str(nrnuc) + "\n")

print(len(processes))
processes = []
logfilearr = []
print(len(processes))

        
for fhccsn in fhccsnarr:
    for fhagb in fhagbarr:
        for inflow in inflowarr:
            filenamenuc = f"Nuclear_fhccsn{fhccsn:.3f}_fhagb{fhagb:.3f}"
            if inflow == True:
                filenamenuc = f"Nuclear_fhccsn{fhccsn:.3f}_fhagb{fhagb:.3f}_inflow"
            filenameglob = f"Global_fhccsn{fhccsn:.3f}_fhagb{fhagb:.3f}"
            with open ("config/"+ suitename +"/" + filenamenuc+".config", "w") as outfile:
                with open ("config/"+ suitename +"/stdnuc_paramtest.config") as infile:
                    outfile.write("output "+outputdir+ "Output/"+suitename + "/" +filenamenuc+"\n")
                    outfile.write("readin-dir "+outputdir+ "/Output/"+suitename + "/"  +filenameglob + "\n")
                    outfile.write(infile.read())
                    outfile.write("\n")
                    outfile.write("fh-ccsn {:.4f}\n".format(fhccsn))
                    outfile.write("fh-nsm {:.4f}\n".format(defaultarr[1]))
                    outfile.write("fh-agb {:.4f}\n".format(fhagb))
                    outfile.write("fh-sn1a {:.4f}\n".format(defaultarr[3]))
                    outfile.write("eject {:.4f}\n".format(defaultarr[5]))
                    outfile.write("inflow-on {}\n".format(inflow))
            logfile = outputdir+ "Output/"+suitename + "/" +filenamenuc+ "/output.log"

            outputfolder = outputdir + "/" + "Output/"+suitename+ "/" + filenamenuc
            if not os.path.exists(outputfolder):
                os.mkdir(outputfolder)
            launchglob = subprocess.Popen("./Ramices_Launch.sh -config config/"+suitename + "/" + filenamenuc + ".config", shell=True, stdout= subprocess.PIPE, bufsize=1)
            processes.append(launchglob)
            logfilearr.append(logfile)

for  nrnuc, (launchnuc, logfile) in enumerate(zip(processes, logfilearr)):
    with open(logfile, 'w') as log:
        for line in iter(launchnuc.stdout.readline, b''):
            #print ("glob" + str(nrglob) + " " + str(line))
            log.write(str(line))
    launchnuc.stdout.close()
    launchnuc.wait()
    print("Done Nuclear Nr " + str(nrnuc) + "\n")

print(len(processes))
processes = []
logfilearr = []
print(len(processes))


for fhccsn in fhccsnarr:
    for fhsn1a in fhsn1aarr:
        for inflow in inflowarr:
            filenamenuc = f"Nuclear_fhccsn{fhccsn:.3f}_fhsn1a{fhsn1a:.4f}"
            if inflow == True:
                filenamenuc = f"Nuclear_fhccsn{fhccsn:.3f}_fhsn1a{fhsn1a:.4f}_inflow"
            filenameglob = f"Global_fhccsn{fhccsn:.3f}_fhsn1a{fhsn1a:.4f}"
            with open ("config/"+ suitename +"/" + filenamenuc+".config", "w") as outfile:
                with open ("config/"+ suitename +"/stdnuc_paramtest.config") as infile:
                    outfile.write("output "+outputdir+ "Output/"+suitename + "/" +filenamenuc+"\n")
                    outfile.write("readin-dir "+outputdir+ "/Output/"+suitename + "/"  +filenameglob + "\n")
                    outfile.write(infile.read())
                    outfile.write("\n")
                    outfile.write("fh-ccsn {:.4f}\n".format(fhccsn))
                    outfile.write("fh-nsm {:.4f}\n".format(defaultarr[1]))
                    outfile.write("fh-agb {:.4f}\n".format(defaultarr[2]))
                    outfile.write("fh-sn1a {:.4f}\n".format(fhsn1a))
                    outfile.write("eject {:.4f}\n".format(defaultarr[5]))
                    outfile.write("inflow-on {}\n".format(inflow))
            logfile = outputdir+ "Output/"+suitename + "/" +filenamenuc+ "/output.log"

            outputfolder = outputdir + "/" + "Output/"+suitename+ "/" + filenamenuc
            if not os.path.exists(outputfolder):
                os.mkdir(outputfolder)
            launchglob = subprocess.Popen("./Ramices_Launch.sh -config config/"+suitename + "/" + filenamenuc + ".config", shell=True, stdout= subprocess.PIPE, bufsize=1)
            processes.append(launchglob)
            logfilearr.append(logfile)

for  nrnuc, (launchnuc, logfile) in enumerate(zip(processes, logfilearr)):
    with open(logfile, 'w') as log:
        for line in iter(launchnuc.stdout.readline, b''):
            #print ("glob" + str(nrglob) + " " + str(line))
            log.write(str(line))
    launchnuc.stdout.close()
    launchnuc.wait()
    print("Done Nuclear Nr " + str(nrnuc) + "\n")

print(len(processes))
processes = []
logfilearr = []
print(len(processes))


for fhccsn in fhccsnarr:
    for ejectglob in ejectglobalarr:
        for inflow in inflowarr:
            filenamenuc = f"Nuclear_fhccsn{fhccsn:.3f}_ejectglob{ejectglob:.3f}"
            if inflow == True:
                filenamenuc = f"Nuclear_fhccsn{fhccsn:.3f}_ejectglob{ejectglob:.3f}_inflow"
            filenameglob = f"Global_fhccsn{fhccsn:.3f}_ejectglob{ejectglob:.3f}"
            with open ("config/"+ suitename +"/" + filenamenuc+".config", "w") as outfile:
                with open ("config/"+ suitename +"/stdnuc_paramtest.config") as infile:
                    outfile.write("output "+outputdir+ "Output/"+suitename + "/" +filenamenuc+"\n")
                    outfile.write("readin-dir "+outputdir+ "/Output/"+suitename + "/"  +filenameglob + "\n")
                    outfile.write(infile.read())
                    outfile.write("\n")
                    outfile.write("fh-ccsn {:.4f}\n".format(fhccsn))
                    outfile.write("fh-nsm {:.4f}\n".format(defaultarr[1]))
                    outfile.write("fh-agb {:.4f}\n".format(defaultarr[2]))
                    outfile.write("fh-sn1a {:.4f}\n".format(defaultarr[3]))
                    outfile.write("eject {:.4f}\n".format(defaultarr[5]))
                    outfile.write("inflow-on {}\n".format(inflow))
            logfile = outputdir+ "Output/"+suitename + "/" +filenamenuc+ "/output.log"

            outputfolder = outputdir + "/" + "Output/"+suitename+ "/" + filenamenuc
            if not os.path.exists(outputfolder):
                os.mkdir(outputfolder)
            launchglob = subprocess.Popen("./Ramices_Launch.sh -config config/"+suitename + "/" + filenamenuc + ".config", shell=True, stdout= subprocess.PIPE, bufsize=1)
            processes.append(launchglob)
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

