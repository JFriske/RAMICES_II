
import subprocess
from datetime import datetime


beg = datetime.now()
#readout time at beginning and output at the end


cgmstartarr = [100, 300, 700]

cgm_pollutionarr = [0.0, 0.5, 1.0]

cgm_coolingarr = [1.0, 10.0, 100.0]

cgmendarr = [7000, 14000]

suitename = "NucDiskCGMTest"

outputdir = "/disk/xray8/jksf/ChemicalEvolution/"
#launchmake = subprocess.Popen("../../make", shell=True, stdout=subprocess.PIPE)
#launchmake.wait()

nrglob = 0
nrnuc = 0 
processes = []
logfilearr = []
for cgm_pollution in cgm_pollutionarr:
    for cgmstart in cgmstartarr:
        for cgm_cooling in cgm_coolingarr:
            
            cgmendmass = 700

            filenameglob = f"Global_cgmstart{cgmstart:.0f}_cgmend{cgmendmass:.0f}_cgmpolluting{cgm_pollution:.2f}_cgmcooling{cgm_cooling:.2f}"
            with open ("config/"+ suitename +"/" + filenameglob+".config", "w") as outfile:
                with open ("config/"+ suitename +"/stdglob_cgmtest.config") as infile:
                    outfile.write("output "+outputdir+ "Output/"+suitename + "/" +filenameglob+"\n")
                    outfile.write(infile.read())
                    outfile.write("\n")
                    outfile.write("cgm-mass {:.4f}\n".format(cgmstart))
                    outfile.write("cgm-mass-end {:.4f}\n".format(700))
                    outfile.write("cgm-pollution {:.4f}\n".format(cgm_pollution))
                    outfile.write("cgm-cooling-factor {:.4f}\n".format(cgm_cooling))

            logfile = outputdir+ "Output/"+suitename + "/" +filenameglob+ "/output.log"

            #nrglob +=1 
            launchglob = subprocess.Popen("./Ramices_Launch.sh -config config/"+suitename + "/" + filenameglob + ".config", shell=True, stdout= subprocess.PIPE, bufsize=1)
            processes.append(launchglob)
            logfilearr.append(logfile)


for cgmendmass in cgmendarr:
    for cgm_pollution in cgm_pollutionarr:
        for cgm_cooling in cgm_coolingarr:
        
            cgmstart = cgmendmass

            filenameglob = f"Global_cgmstart{cgmstart:.0f}_cgmend{cgmendmass:.0f}_cgmpolluting{cgm_pollution:.2f}_cgmcooling{cgm_cooling:.2f}"
            with open ("config/"+ suitename +"/" + filenameglob+".config", "w") as outfile:
                with open ("config/"+ suitename +"/stdglob_cgmtest.config") as infile:
                    outfile.write("output "+outputdir+ "Output/"+suitename + "/" +filenameglob+"\n")
                    outfile.write(infile.read())
                    outfile.write("\n")
                    outfile.write("cgm-mass {:.4f}\n".format(cgmstart))
                    outfile.write("cgm-mass-end {:.4f}\n".format(cgmendmass))
                    outfile.write("cgm-pollution {:.4f}\n".format(cgm_pollution))
                    outfile.write("cgm-cooling-factor {:.4f}\n".format(cgm_cooling))


            logfile = outputdir+ "Output/"+suitename + "/" +filenameglob+ "/output.log"

            launchglob = subprocess.Popen("./Ramices_Launch.sh -config config/"+suitename + "/" + filenameglob + ".config", shell=True, stdout=  subprocess.PIPE, bufsize=1)
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

for cgm_pollution in cgm_pollutionarr:
    for cgmstart in cgmstartarr:
        for cgm_cooling in cgm_coolingarr:

            cgmendmass =  700
            
            filenamenuc = f"Nuclear_cgmstart{cgmstart:.0f}_cgmend{cgmendmass:.0f}_cgmpolluting{cgm_pollution:.2f}_cgmcooling{cgm_cooling:.2f}"
            with open ("config/"+ suitename +"/" +filenamenuc+ ".config", "w") as outfile:
                with open("config/"+ suitename +"/stdnuc_cgmtest.config") as infile:
                    outfile.write("output "+outputdir+ "Output/"+suitename + "/" +filenamenuc+"\n")
                    outfile.write("readin-dir "+outputdir+ "Output/"+suitename + "/"  +filenameglob + "\n")
                    outfile.write(infile.read())

            logfile = outputdir+ "Output/"+suitename + "/" +filenamenuc+ "/output.log"
            logfilearr.append(logfile)

            launchnuc = subprocess.Popen("./Ramices_Launch.sh -config config/"+suitename + "/" + filenamenuc + ".config", shell=True, stdout= subprocess.PIPE, bufsize=1)
            processes.append(launchnuc)

for cgmendmass in cgmendarr:
    for cgm_pollution in cgm_pollutionarr:
        for cgm_cooling in cgm_coolingarr:
        
            cgmstart = cgmendmass
            filenamenuc = f"Nuclear_cgmstart{cgmstart:.0f}_cgmend{cgmendmass:.0f}_cgmpolluting{cgm_pollution:.2f}_cgmcooling{cgm_cooling:.2f}"
            with open ("config/"+ suitename +"/" +filenamenuc+ ".config", "w") as outfile:
                with open("config/"+ suitename +"/stdnuc_cgmtest.config") as infile:
                    outfile.write("output "+outputdir+ "Output/"+suitename + "/" +filenamenuc+"\n")
                    outfile.write("readin-dir "+outputdir+ "Output/"+suitename + "/"  +filenameglob + "\n")
                    outfile.write(infile.read())
            
            logfile = outputdir+ "Output/"+suitename + "/" +filenamenuc+ "/output.log"
            logfilearr.append(logfile)
            
            launchnuc = subprocess.Popen("./Ramices_Launch.sh -config config/"+suitename + "/" + filenamenuc + ".config", shell=True, stdout= subprocess.PIPE, bufsize=1)
            processes.append(launchnuc)


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
