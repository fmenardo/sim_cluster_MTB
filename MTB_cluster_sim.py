#sim_cluster_MTB a pipeline to simulate the molecular evolution of MTB epidemics, compute clustering rates and terminal branch lengths
#Copyright (C) 2022  Fabrizio Menardo

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <https://www.gnu.org/licenses/>.




from Bio import SeqIO
from ete3 import Tree
import argparse
import os
import csv
import sys
import re
import subprocess
import random
import shutil


BEAST_PATH ="~/software/beast/bin/beast"
#SCRATCH_PATH=  "path/to/scratch"
SCRATCH_PATH= os.getcwd()


def generate_xml(lin,br,er,dr,sr,stem,sim_time,end_condition,min_mt,max_mt):

    xml = """<beast version='2.0' namespace='master:master.model:master.steppers:master.conditions:master.postprocessors:master.outputs'>
  <run spec='InheritanceTrajectory' samplePopulationSizes='False'"""+str(sim_time)+""">

    <model spec='Model' id='model'>
      <population spec='Population' id='E' populationName='E'/>
      <population spec='Population' id='I' populationName='I'/>
      <population spec='Population' id='S' populationName='S'/>

      <reaction spec='Reaction' reactionName='infection' rate='"""+str(br)+"""'>
	 I:1 -> I:1 + E:1
      </reaction>

      <reaction spec='Reaction' reactionName='becom_infectious' rate='"""+str(er)+"""'>
	 E:1 -> I:1 
      </reaction>

      <reaction spec='Reaction' reactionName='death' rate='"""+str(dr)+"""'>
	I -> 0
      </reaction>

      <reaction spec='Reaction' reactionName='sample' rate='"""+str(sr)+"""'>
	I:1 -> S:1
      </reaction>

    </model>
    
    <initialState spec='InitState'>
      <lineageSeed spec='Individual' population='@I'/>

    </initialState>

    """+str(end_condition)+"""
 
     <inheritancePostProcessor spec='LineageFilter'
        populationName='E'
        discard='true'
        noClean='false'
        leavesOnly='false' />


     <inheritancePostProcessor spec='LineageFilter'
        populationName='I'
        discard='true'
        noClean='false'
        leavesOnly='false' />

    <inheritancePostProcessor spec='LineageFilter'
        populationName='S'
        discard='true'
        noClean='false'
        leavesOnly='true' />

    <postSimCondition spec='LeafCountPostSimCondition' 
                      nLeaves='"""+str(min_mt)+"""'
                      exact ='false'
                      exceedCondition='true'/>
   
    <postSimCondition spec='LeafCountPostSimCondition' 
                      nLeaves='"""+str(max_mt)+"""'
                      exact ='false'
                      exceedCondition='false'/>
    
    <output spec='NewickOutput' fileName='"""+stem+""".newick' collapseSingleChildNodes='true'/>
    <output spec='NewickOutput' fileName='"""+stem+""".tr.newick'/>
  </run>
</beast>"""

    F=open(stem+".xml","w")
    F.write(xml)
    F.close()
    return()


def prune_t(leaf_to_prune,tree):

	G = tree.search_nodes(name=leaf_to_prune)[0]		
	parent= G.up
	dist_parent=G.dist

	if (len(parent.get_children()) == 2):
	

		if  parent.children[0] != G:
			parent.children[0].dist = parent.children[0].dist + parent.dist		

		if parent.children[1] != G:
			parent.children[1].dist = parent.children[1].dist + parent.dist	
			
	G.detach()

	if (len(parent.get_children()) == 1):
		parent.delete()		# after pruning the remaining branch will be like this ---/---leaf_name. I delete useless node keeping the b length


	return (tree)

def parse_tree (time_sampling, stem):
    leaves_dict={}
    t_height=0

    t = Tree(stem + ".newick",format=1)

    ori_t = t.copy()
    for node in ori_t.traverse():
    	if node.is_root():
    		root=node
    		root.name="root"
    		MR_leaf=root.get_farthest_leaf()

    		t_height=round(MR_leaf[1],10)

    min_height = t_height - time_sampling

    leaves=t.get_leaves()

    for leaf in (leaves):
    	l_height= round(ori_t.get_distance("root",leaf.name),10)

    	if l_height < min_height:
        	t = prune_t(leaf.name,t)

###  add outgroup and fix weird node with single child 
    t1=Tree()
    t1.add_child(name="ANC",dist="1")
    t1.add_child(t)
    t=t1

    for node in t.traverse():
    	if node.is_root():
    		root=node
    		root.name="root1"
    		MR_leaf=root.get_farthest_leaf()
    		t_height=round(MR_leaf[1],10)

    	
    for ch in range(0,2):
    	if  root.children[ch].name != "ANC" and len(root.children[ch].get_children()) == 1:
    		child=root.children[ch]
    		child.children[0].dist = child.children[0].dist + child.dist
    		child.delete()	


    F=open(stem+".newick_sampled","w")
    F.write(t.write(format=5))
    F.close()

    leaves=t.get_leaves()

    for leaf in (leaves):
    	l_height= round(t.get_distance("root1",leaf.name),10)
    	leaves_dict.update({leaf.name : t_height-l_height})

    with open(stem+'_lname_time.csv', 'w') as f:  
        writer = csv.writer(f)
        for k, v in leaves_dict.items():
            writer.writerow([k, v])

    return(leaves_dict,t_height)


def rescale_tree(stem_cr,invariant):

	var=4000000-invariant
	print("rescaling tree")

	t = Tree(stem_cr + ".fasta_var.raxml.bestTree",format=1)

	for node in t.traverse():

		if (node.dist < 0.00000001):
			node.dist=0
		else:
			node.dist= (node.dist*var)/4000000


	F=open(stem_cr + ".raxml.rescaled","w")
	F.write(t.write(format=5))
	F.close()

	for node in t.traverse():

		node.dist= round(node.dist*4000000)

	prune_t("ANC",t)

	F=open(stem_cr + ".raxml.rescaled_SNP","w")
	F.write(t.write(format=5))
	F.close()



	return(t)

def run_snp_site (stem_cr, folder_cr, cr):

        print("running snp-site clock rate " + str(cr))
        
        command= "snp-sites -o " + stem_cr + ".fasta_var " +  folder_cr + "/" + stem_cr + ".fasta"
        p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, error = p.communicate()

        #out=os.system("snp-sites -o " + stem_cr + ".fasta_var " +  folder_cr + "/" + stem_cr + ".fasta")
        snpsite_err=b"Warning: No SNPs were detected"
        if error.find(snpsite_err) >=0 :
            sys.exit("this combination of parameters resulted in no SNPs, consider increasing the number of strains to sample, or the molecular clock rate")


        os.system("snp-sites -C -o " + stem_cr + ".fasta_count " +  folder_cr + "/" + stem_cr + ".fasta")

        with open(stem_cr+'.fasta_count', 'r') as inv_sites:
                reader = csv.reader(inv_sites)
                for row in reader:
                        A=row[0]
                        C=row[1]
                        G=row[2]
                        T=row[3]
                inv_sites.close()

        return(A,C,G,T)

def clustering(t,stem_cr,SNP_t):

        print ("clustering")
        perc=[]
        normalized_cluster_n=[]
        leaves=t.get_leaves()

        for x in range(0, SNP_t+1):
                clus=[]
                os.system("TreeCluster.py -t " + str(x) + " -m max -i " + stem_cr + ".raxml.rescaled_SNP -o " + stem_cr + ".cl" + str(x))

                file = open(stem_cr + ".cl" + str(x) ,"r")


                data = file.read()

                occurrences = data.count("-1")
               
                prop= round(1-(int(occurrences) / int(len(leaves))),4)
                perc.append(str(prop))
                file.close()
                #print (perc[x])

                file = open(stem_cr + ".cl" + str(x) ,"r")
                Lines = file.readlines()

                for line in Lines:

                    cluster = re.match('^\d+\s(\d+)', line)

                    if (cluster):
                        clus.append(cluster.group(1))
                clus = list(dict.fromkeys(clus))

                n_clus=round(int(len(clus))/int(len(leaves)),4)

                normalized_cluster_n.append(str(n_clus))
                #print (str(n_clus))
                file.close()
        perc.append("\n")
        F=open(stem_cr + ".cl_r","w")
        F.write("	".join(perc))
        F.close()

#        normalized_cluster_n.append("\n")
#        F=open(stem_cr + ".clnn_r","w")
#        F.write("	".join(normalized_cluster_n))
#        F.close()


        return(leaves)

def terminal_b_length(leaves,stem_cr):

        leaves_dist= {}

        for leaf in leaves :

                leaf.dist
                leaves_dist.update({leaf.name : leaf.dist})

        with open(stem_cr +'_ldist.csv', 'w') as f:
                writer = csv.writer(f)
                for k, v in leaves_dist.items():
                        writer.writerow([k, v])


def subsample_fasta_var (s,stem_cr,stem_cr1,t_height,ps_sr):


	tot={}
	length=[[]for _ in range(len(ps_sr))]
	leaves_temp=[[]for _ in range(len(ps_sr))]
	OUT=[{}for _ in range(len(ps_sr))]
	out_l=[{}for _ in range(len(ps_sr))]
	
	invA=[]
	invT=[]
	invC=[]
	invG=[]

	years=s.split(",")


	for k, v in leaves_dict.items():

		year_flag=0
		for year in years:
			if  ((int(year)-1) <= v) and (v < int(year)):
				year_flag=1
		if year_flag ==1:

			for x in range(0,len(ps_sr)):

				if random.random() <= float(ps_sr[x]):
					
					leaves_temp[x].append(k)


	for x in range(0, len(ps_sr)):
		for record in SeqIO.parse(stem_cr1 + ".fasta_var", "fasta"):
			length[x] = len(record.seq)
			flag=0
			for leaf in leaves_temp[x]:
				if (leaf == str(record.id)):
					flag=1
				if ("ANC" == str(record.id)):
					flag=1
			if flag==1:
				out_l[x].update({str(record.id):str(record.seq)})
				OUT[x].update({str(record.id):""})

	for x in range(0, len(ps_sr)):


		N=0
		count_variable=0
		count_invariant=0
		count_missing=0
		count_invA=0
		count_invT=0
		count_invC=0
		count_invG=0
		list_pos=[]

		while N < length[0]:
	
			(OUT[x],message)=find_variable(N,out_l[x],OUT[x])

			N=N+1

			if message == "variable":
				count_variable =int(count_variable)+1
				list_pos.append(N)

			if message == "missing":
				count_missing=int(count_missing) +1	
			if message == "A":	
				count_invariant = int(count_invariant)+1
				count_invA = int(count_invA) + 1	
			if message == "T":	
				count_invariant = int(count_invariant)+1
				count_invT = int(count_invT) + 1
			if message == "C":	
				count_invariant = int(count_invariant)+1
				count_invC = int(count_invC) + 1
			if message == "G":	
				count_invariant = int(count_invariant)+1
				count_invG = int(count_invG) + 1

		file_out = stem_cr + "_" +str(ps_sr[x]) + "_" + str(sim_id) + ".fasta_var"

		with open(file_out, 'w') as f_out:
			for seq_name, sequence in OUT[x].items():
				f_out.write(">"+ str(seq_name) +"\n" + str(sequence) + "\n")
		f_out.close()

		invA.append(count_invA)
		invT.append(count_invT)
		invC.append(count_invC)
		invG.append(count_invG)
	
	return(invA, invC, invG, invT)

	
def find_variable(N,out_l,OUT):
	
	message=""
	
	pos=""
	for name, seq in out_l.items():
		pos += seq[N]
	
	AA=len(re.findall('[Aa]', pos))		
	TT=len(re.findall('[Tt]', pos))
	GG=len(re.findall('[Gg]', pos))
	CC=len(re.findall('[Cc]', pos))
	if (AA > 0):
		AA=1
	if (TT > 0):
		TT=1
	if (CC > 0):
		CC=1
	if (GG > 0):
		GG=1
	if ((AA+CC+GG+TT) > 1):
		message = "variable"
			
		#print "due"
		#print A + T + G +C
		for name_out, seq_out in out_l.items():
			OUT[name_out] += seq_out[N]
	else:
		if AA == 1:
			message="A"
		if TT == 1:
			message="T"
		if CC == 1:
			message="C"
		if GG == 1:
			message="G"

	return (OUT,message)

def run_raxml(stem_cr):

	command= "raxml-ng --redo --thread 1 --blmin 0.0000000001 --outgroup ANC --model HKY --tree rand{" + str(rrt) + "},pars{" + str(rpt) + "} --msa " +  stem_cr + ".fasta_var"
	p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	out, error = p.communicate()


	raxml_err=b"ERROR: Your alignment contains less than 4 sequences"
	if out.find(raxml_err) >=0 :
		sys.exit("this combination of parameters resulted in less than 4 taxa in the raxml analysis, consider increasing the number of strains to sample")


	raxml_err=b"Alignment sites / patterns: 1 / 1"
	if out.find(raxml_err) >=0 :
		sys.exit("this combination of parameters resulted 1 single SNP in the raxml analysis, consider increasing the number of strains to sample or the clock rate")

	raxml_err=b"ERROR: Error loading MSA: cannot parse any format supported by"
	if out.find(raxml_err) >=0 :
		sys.exit("this combination of parameters resulted no SNPs in the raxml analysis, consider increasing the number of strains to sample or the clock rate")

	raxml_err=b"Alignment comprises 0 taxa, 0 partitions and 0 patterns"
	if out.find(raxml_err) >=0 :
		sys.exit("alignment error (0 taxa, 0 partitions and 0 patterns)  in the raxml analysis, consider increasing the number of strains to sample")

	raxml_err=b"Alignment sites / patterns: 2 /"
	if out.find(raxml_err) >=0 :
		sys.exit("this combination of parameters resulted in only 2 SNPs in the raxml analysis, consider increasing the number of strains to sample")

	raxml_err=b"Alignment sites / patterns: 3 /"
	if out.find(raxml_err) >=0 :
		sys.exit("this combination of parameters resulted in only 3  SNPs in the raxml analysis, consider increasing the number of strains to sample")
	raxml_err=b"Alignment sites / patterns: 4 /"
	if out.find(raxml_err) >=0 :
		sys.exit("this combination of parameters resulted in only 4  SNPs in the raxml analysis, consider increasing the number of strains to sample")


			#os.system("raxml-ng --threads 1 --redo --blmin 0.0000000001 --outgroup ANC --model HKY --msa " +  stem_cr + ".fasta_var")
# wrapper

parser = argparse.ArgumentParser()

#parser.add_argument('INFILE',type=str,help='path to the newick tree')
parser.add_argument('-l','--lineages', metavar='lineages', default='10', help='minimum number of living lineages (I population), when the simulation exceed this number it stops' , type =int)
parser.add_argument('-ts','--time_sampling', metavar='time', default='', help='number of years of sampling (starting from present and going backward. default=(10)' , type =int, nargs=1)
parser.add_argument('-br','--birth_rate', metavar='B_R', default='', help='birth rate' , type =float, nargs=1)
parser.add_argument('-dr','--death_rate', metavar='D_R', default='', help='death rate (excluding sampling)' , type =float, nargs=1)
parser.add_argument('-sr','--sampling_rate', metavar='S_R', default='', help='sampling rate (excluding death, if death rate = 0, this is the death rate, otherwise the death rate is S_R + D_R' , type =float, nargs=1)
parser.add_argument('-er','--exposed_rate', metavar='E_R', default='', help='rate at which exposed become infectious' , type =float, nargs=1)
parser.add_argument('-sim_n','--simulation_number', metavar='SIM', default='', help='simulation number ID' , type =int, nargs=1)
parser.add_argument('-cr','--clock_rate', metavar='C_R', default='', help='clock rate  (nucleotide substitution per site per year), multiple values possible' , type =float, nargs='+')
parser.add_argument('-ps_sr','--post_sim_sampling_rates', metavar='PS_SR', default='', help='probability of each strain to be sampled (post simulation), multiple rates possible at once: eg. <-ps_sr 1 0.5 0.1', type = str, nargs='*')
parser.add_argument('-ps_sy','--post_sim_sampling_years', metavar='PS_SY', default='', help='sample only in these years (post sim), multiple scheme possible at once: eg. <-ps_sy 1,2,3 1,3,5>  default (all) ', type = str, nargs='*',)
parser.add_argument('-c','--clean', default=False, help='delete all intermediate file, keep only clustering results and terminal branch lengths (default: False)',action='store_true')
parser.add_argument('-s','--stop', default=False,metavar="lineages|time", help='stop criterion, the MASTER simulation should stop when reaching a certain number of infectious existing lineages("lineages"; specified with -l) or after a certain time ("time", specified with -t)(default = "lineages"',type=str,choices=["lineages","time"])
parser.add_argument('-min_mt','--min_master_tips', metavar='', default='4', help='minimum number of tips in the tree output of MASTER' , type =int)
parser.add_argument('-max_mt','--max_master_tips', metavar='', default='2500', help='max number of tips in the tree output of MASTER' , type =int)
parser.add_argument('-t','--time', metavar='', default='10', help='Simulation time for MASTER, to be used with "time" as stop criterion ' , type =int)
parser.add_argument('-rpt','--raxml_pars_tree', metavar='', default='1', help='number of starting parsimony trees for raxml' , type =int)
parser.add_argument('-rrt','--raxml_rand_tree', metavar='', default='1', help='number of starting random trees for raxml' , type =int)
parser.add_argument('-SNP_t','--SNP_threshold', metavar='', default='50', help='clustering will be performed for all values in the interval 0-SNP_t' , type =int)
parser.add_argument('-f','--force', metavar='', default='0', help='Discard simulation if tree height < f (default: 0)',type =int)



arguments = parser.parse_args()

lin=arguments.lineages
time_sampling=arguments.time_sampling[0]
br=arguments.birth_rate[0]
dr=arguments.death_rate[0]
sr=arguments.sampling_rate[0]
er=arguments.exposed_rate[0]
sim_id=arguments.simulation_number[0]
c_rates=arguments.clock_rate
ps_sy=arguments.post_sim_sampling_years
ps_sr=arguments.post_sim_sampling_rates
min_mt=arguments.min_master_tips
max_mt=arguments.max_master_tips
time=arguments.time
rpt=arguments.raxml_pars_tree
rrt=arguments.raxml_rand_tree
SNP_t=arguments.SNP_threshold
force_t=arguments.force


if (arguments.stop == False):
	stop="lineages"
else:
	stop = arguments.stop

print(stop)

if (stop == "time"):

	sim_time = (" simulationTime=\"" + str(time)+"\"")
	lin = 0
	end_condition=""
if (stop == "lineages"):

	sim_time =""
	time = 0
	end_condition="""<lineageEndCondition spec='LineageEndCondition' 
                        nLineages='0'
                        isRejection='true'
                        population="@I"/>

    <lineageEndCondition spec='LineageEndCondition' 
                        nLineages='"""+str(lin)+"""'
                        isRejection='false'
                        alsoGreaterThan="True"
                        population="@I"/>"""


mt=(str(min_mt)+"-"+str(max_mt))

stem= str(br) + "_" + str(er) + "_" + str(dr) + "_" + str(sr)  + "_" + str(lin) + "_" + str(mt) + "_" + str(time_sampling) + "_" + str(time) + "_" + str(force_t) + "_" + str(sim_id)
folder = "sim_"+  str(br) + "_" + str(er) + "_" + str(dr) + "_" + str(sr)  + "_" + str(lin) + "_" + str(mt) + "_" + str(time_sampling) + "_" + str(time) + "_" + str(force_t)


# create folder 
if not os.path.exists(folder):
	os.makedirs(folder)
os.chdir(folder)

# generate xml for MASTER

generate_xml(lin,br,er,dr,sr,stem,sim_time,end_condition,min_mt,max_mt)


flag_force=0
while (flag_force == 0):

# run MASTER

    os.system(BEAST_PATH + " -noerr " + stem + ".xml")  #-noerr 

#parse tree

    leaves_dict, t_height = parse_tree (time_sampling, stem)

#force t height

    if(t_height < force_t):
        print ("tree height < than sampling time (" + str(t_height) + " < " + str(force_t) + "rejecting_simulation")
        continue
        

## for each clock rate

    for n in range(0, len(c_rates)):
        cr=c_rates[n]



# create folder in scratch

        stem_cr= str(br) + "_" + str(er) + "_" + str(dr) + "_" + str(sr)  + "_" + str(lin) + "_" + str(mt) + "_" + str(time_sampling) + "_" + str(time) + "_" + str(force_t) + "_" + str(cr) + "_" + str(sim_id)
        folder_cr= SCRATCH_PATH + "/" + folder + "/" + str(br) + "_" + str(er) + "_" + str(dr) + "_" + str(sr)  + "_" + str(lin) + "_" + str(mt) + "_" + str(time_sampling) + "_" + str(time) + "_" + str(force_t) + "_" + str(cr)

        if not os.path.exists(folder_cr):
        	os.makedirs(folder_cr)


# run seq-gen

        os.system("seq-gen -a 1 -l 4000000 -m HKY -t 2 -of -s " + str(cr) + " "  + stem + ".newick_sampled > " + folder_cr + "/" + stem_cr + ".fasta")

# create folder in wd

        folder_cr_wd=  str(br) + "_" + str(er) + "_" + str(dr) + "_" + str(sr)  + "_" + str(lin) + "_" + str(mt) + "_" + str(time_sampling) + "_" + str(time) + "_" + str(force_t) + "_" + str(cr)

        if not os.path.exists(folder_cr_wd):
                    os.makedirs(folder_cr_wd)

        os.chdir(folder_cr_wd)

# keep only variable sites, but count nonvariable ones for asc correction

        (A_o,C_o,G_o,T_o) = run_snp_site(stem_cr, folder_cr, cr)

# subsampling only the strains from the specified years and create additional fasta_var files, also run variable fasta and integarte inv count with the count of snp sites

        flag_subsampling = ""

        if ps_sy == "":
            for X in range(1,(int(time_sampling)+1)):
                if X ==1:
                    string = str(X)
                else:
                    string= string + "," + str(X)
            ps_sy =[string]
            flag_subsampling = 'PASS'

        stem_cr1 =stem_cr

        if ps_sr == "":
            ps_sr =[1]



        for s in ps_sy:
            A = T = C = G = 0



            stem_cr= str(br) + "_" + str(er) + "_" + str(dr) + "_" + str(sr)  + "_" + str(lin) + "_" + str(mt) + "_" + str(time_sampling) + "_" + str(time) + "_" + str(force_t) + "_" + str(cr) + "_" + str(s) 
            (invA, invC, invG, invT) = subsample_fasta_var(s,stem_cr,stem_cr1,t_height,ps_sr)

            ktr=0
            stem_cr2=stem_cr

            for pssr in ps_sr:

                stem_cr= str(br) + "_" + str(er) + "_" + str(dr) + "_" + str(sr)  + "_" + str(lin) + "_" + str(mt) + "_" + str(time_sampling) + "_" + str(time) + "_"+ str(force_t) + "_" + str(cr) + "_" + str(s) + "_" + str(pssr) + "_" + str(sim_id)
                A = int(A_o) + int(invA[ktr])
                C = int(C_o) + int(invC[ktr])
                G = int(G_o) + int(invG[ktr])
                T = int(T_o) + int(invT[ktr])
                invariant = A+C+T+G
                variant = 4000000-invariant
                file_SNP = open(stem_cr + ".SNP_count", "w")
                file_SNP.write(str(variant))
                file_SNP.close()

# make tree

                print ("running raxml clock rate " +str(cr))

                run_raxml(stem_cr)


# rescale tree to have branch lengths in SNPs and remove ANC

                t = rescale_tree(stem_cr, invariant)

# cluster

                leaves = clustering(t,stem_cr,SNP_t)


# calculate lengths terminal branches

                terminal_b_length(leaves,stem_cr)

# mv results files

                os.system("mv " + stem_cr + ".cl_r ../.")
                #os.system("mv " + stem_cr + ".clnn_r ../.")
                os.system("mv " + stem_cr + "_ldist.csv ../.")
                os.system("mv " + stem_cr + ".fasta_var ../.")
                #os.system("mv " + stem_cr + ".SNP_count ../.")
                os.system("mv " + stem_cr + ".raxml.rescaled ../.")

                ktr=ktr+1
# change directory ../. for next clock rate

        os.chdir("../.")


# clean files

        if arguments.clean:
            shutil.rmtree(folder_cr)
            if os.path.exists(folder_cr_wd):
                shutil.rmtree(folder_cr_wd)


    if arguments.clean:
        os.remove(stem + "_lname_time.csv")
        os.remove(stem + ".newick_sampled")
        os.remove(stem + ".xml")
        new=folder_cr.rsplit('/', 1)
        if os.getcwd() != new[0]:
            shutil.rmtree(new[0])

    flag_force=1
