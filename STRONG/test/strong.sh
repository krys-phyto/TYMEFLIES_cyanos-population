################################
### Krys Kibler
### Purpose: Use strong to delineate strains
### Resource: https://gitpub.wei.wisc.edu/htcondor/strong-strain-resolution-docker


################################
### Submitting to the condor job cluster
condor_submit STRONG_submit.submit

# need to edit condor job after submission or it will stall at idle
condor_qedit 1423031 requirements='(TARGET.COLLECTOR_HOST_STRING == "scarcity-cm.glbrc.org:9618?sock=collector") && (Target.OpSysandVer == "CentOS7") && TARGET.HasDocker && (TARGET.Disk >= RequestDisk) && (TARGET.Memory >= RequestMemory) && (TARGET.Cpus >= RequestCpus)'


# View progress of the job
condor_tail -stderr 1364911

# View current data output without interrupting job
condor_ssh_to_job 1364911
  # type "exit" when satisfied


################################
### Batching sample files for multiple STRONG runs

### Just 2015
"rd0216" "rd0665" "rr0007" "rr0052" "rr0055" "rr0059" "rr0068" "rr0088" "rr0098" "rr0105" "rr0111"
"rr0122" "rr0134" "rr0135" "rr0143" "rr0166" "rr0174" "rr0184" "rr0283" "rr0294" "rr0313" "rr0370"
"rr0390" "rr0395" "rr0402" "rr0416" "rr0428" "rr0442" "rr0468" "rr0469" "rr0470" "rr0489" "rr0493"
"rr0543" "rr0611" "rr0628" "rr0646" "rr0670" "rr0707" "rr0769" "rr0793" "rr0798" "rr0801" "rr0838"
"rr0844" "rr0870" "rr0873" "rr0880" "rr0891" "rr0898" "rr0910" "rr0911" "rr0932" "rr0960" "rr0981"
"rr0991"

# put these names into a txt file and create for loop to move the ones that are
# them into samples folder for strong to run

while read folder
do
      mv metags/$folder samples/
done < batch_list_2015.txt
