
import subprocess
from datetime import datetime


beg = datetime.now()
#readout time at beginning and output at the end

# global changes needed
fhsn1aarr = [0.95, 0.99, 0.999]
fhccsnarr = [0.3, 0.6, 0.9]
fhnsmarr = [0.65, 0.8, 0.95]
fhagbarr = [0.1, 0.35, 0.6]

#only global 
ejectglobalarr = [0.3, 0.45, 0.6]

# only nuclear
ejectnucleararr = [0.45, 0.75, 0.95]
hottransferlossarr = [0.7]
coldtransferlossarr = [0.5]

cgmstartarr = [100, 300, 700]

McMillanInflowarr = [True, False]


suitename = "NucDiskCGMParams"

outputdir = "/disk/xray8/jksf/ChemicalEvolution/"
#launchmake = subprocess.Popen("../../make", shell=True, stdout=subprocess.PIPE)
#launchmake.wait()

nrglob = 0
nrnuc = 0 

for fhsn1a in fhsn1aarr:
    for fhccsn in fhccsnarr:
        for fhnsm in fhnsmarr:
            for fhagb in fhagbarr:
                for ejectglob in ejectglobalarr:
                    for cgmstart in cgmstartarr:
                        for McMillanInflow in McMillanInflowarr:
                            filenameglob = f"Global_fhsn1a{fhsn1a:.3f}_fhccsn{fhccsn:.3f}_fhnsm{fhnsm:.3f}_fhagb{fhagb:.3f}_ejectglob{ejectglob:.3f}_cgmstart{cgmstart:.0f}_inflowMcMillan{McMillanInflow}"
                            with open ("config/"+ suitename +"/" + filenameglob+".config", "w") as outfile:
                                with open ("config/"+ suitename +"/globalbaseconfig.config") as infile:
                                    outfile.write("output "+outputdir+ "Output/"+suitename + "/" +filenameglob+"\n")
                                    outfile.write(infile.read())
                                    outfile.write("fh-sn1a {:.4f}\n".format(fhsn1a))
                                    outfile.write("fh-ccsn {:.4f}\n".format(fhccsn))
                                    outfile.write("fh-nsm {:.4f}\n".format(fhnsm))
                                    outfile.write("fh-agb {:.4f}\n".format(fhagb))
                                    outfile.write("eject {:.4f}\n".format(ejectglob))
                                    outfile.write("cgm-mass {:.4f}\n".format(cgmstart))
                                    if (McMillanInflow):
                                        outfile.write("M0 {:.4f}\n".format(0.1))
                                        outfile.write("M1 {:.4f}\n".format(50))
                                        outfile.write("M2 {:.4f}\n".format(100))
                                        outfile.write("b1 {:.4f}\n".format(1))
                                        outfile.write("b2 {:.4f}\n".format(9))
                                    else:
                                        outfile.write("M0 {:.4f}\n".format(0.1))
                                        outfile.write("M1 {:.4f}\n".format(4))
                                        outfile.write("M2 {:.4f}\n".format(80))
                                        outfile.write("b1 {:.4f}\n".format(0.3))
                                        outfile.write("b2 {:.4f}\n".format(8))

                            nrglob +=1 
                            launchglob = subprocess.Popen("./Ramices_Launch.sh -config config/"+suitename + "/" + filenameglob + ".config", shell=True, stdout=subprocess.PIPE, bufsize=1)
                            for line in iter(launchglob.stdout.readline, b''):
                                print ("glob" + str(nrglob) + " " + str(line))
                            launchglob.stdout.close()
                            launchglob.wait()


                            for ejectnuc in ejectnucleararr:
                                for hottransferloss in hottransferlossarr:
                                    for coldtransferloss in coldtransferlossarr:
                                        filenamenuc = f"Nuclear_fhsn1a{fhsn1a:.3f}_fhccsn{fhccsn:.3f}_fhnsm{fhnsm:.3f}_fhagb{fhagb:.3f}_ejectglob{ejectglob:.3f}_ejectnuc{ejectglob:.3f}_cgmstart{cgmstart:.0f}_inflowMcMillan{McMillanInflow}"
                                        with open ("config/"+ suitename +"/" +filenamenuc+ ".config", "w") as outfile:
                                            with open("config/"+ suitename +"/nucbaseconfig.config") as infile:
                                                outfile.write("output "+outputdir+ "Output/"+suitename + "/" +filenamenuc+"\n")
                                                outfile.write("readin-dir "+outputdir+ "Output/"+suitename + "/"  +filenameglob + "\n")
                                                outfile.write(infile.read())
                                                outfile.write("fh-sn1a {:.4f}\n".format(fhsn1a))
                                                outfile.write("fh-ccsn {:.4f}\n".format(fhccsn))
                                                outfile.write("fh-nsm {:.4f}\n".format(fhnsm))
                                                outfile.write("fh-agb {:.4f}\n".format(fhagb))
                                                outfile.write("eject {:.4f}\n".format(ejectnuc))
                                                outfile.write("hot-gas-transport-loss {:.4f}\n".format(hottransferloss))
                                                outfile.write("cold-gas-transport-loss {:.4f}\n".format(coldtransferloss))
                                        
                                        nrnuc +=1 
                                        launchnuc = subprocess.Popen("./Ramices_Launch.sh -config config/"+suitename + "/" + filenamenuc + ".config", shell=True, stdout=subprocess.PIPE, bufsize=1)
                                        for line in iter(launchnuc.stdout.readline, b''):
                                            print ("nuc" + str(nrnuc) + " " + str(line))
                                        launchnuc.stdout.close()
                                        #print (launchnuc.communicate()[0])
                                        launchnuc.wait()

                                
end = datetime.now() -beg

current_time = end#.strftime("%H:%M:%S")
print("Current Time =", current_time)





# fh-sn1a 0.99
# fh-ccsn 0.7
# fh-nsm 0.4
# fh-agb 0.5


# output Output/NucDiskFindHook

# hot-gas-transport-loss 0.7
# cold-gas-transport-loss 0.5

# eject 0.4 (one for nuc and one for total)

