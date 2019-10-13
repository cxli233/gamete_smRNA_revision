#This uses the information in a BS Seeker 2 ATCGmap file to determine the mean methylation level in and near every feature in a GFF file.
#The exact range to be analyzed needs to be set by the user in variables below
#GFF file must be sorted numerically, chromosome than start position. 
#Overlapping GFFs are tolerated as long as they are sorted.
#GFF file should only include numeric chromosomes, no scaffolds or organellar.
#ATCGmap file must be sorted numerically, chromosome than position (sort -n -k1,1 -k3,3). 
#All non-numeric chromosomes are skipped. 

#Derived from BS2_GFF_range_methylator5.py. 
#The major difference is that CCH is also measured as a subset of total CHH. 

#OUTLINE
#Create a set of bins corresponding to specific distances from each feature.
#Get data from ATCGmap.sorted file for each position that in a bin using Bins_Meth_Analyzer, which contains a nested GFF_line_extractor function.
	#Record the number of cytosines in each bin and the number read as methylated or unmethylated
	#Calculate summary stats for each bin for each methylation context.
	#If there are no overlapping features, the above will occur in a single pass through ATCGmap file for each bin. 
	#If there are overlapping features, multiple runs through will be made the ATCGmap file, each time with a smaller set of features for each bin
    	#(the "residuals" or overlapping feature that were missed in the prior run)
#Repeat above for every input file in directory

################################################################################################################
#INITIALIZE SOME VARIABLES
################################################################################################################


#I/O

#Directory of ATCGmaps (unlike previous versions, don't include '/' at end!)
partial_ATCGmap_path = '/scratch/gent/rice_BS/sorted_maps'

#ATCGmap file endings
ending = "7.ATCGmap.sorted"

#GFF file
GFF_path = '/home/gent/references/MSU7_filtered_genes.gff' #File with coordinates in GFF format

#output ending
output_ending = "_MSU7_filtered_genes_upstream.meths"


#Bin specifications 
upstream_edge = True #Define whether bins are defined by upstream or downstream edges of features	
bin_size = 100 #length of bin in bp
bin_start = -3000 #distance from start position of first bin to edge of feature. "-" means upstream. 
				  #e.g., for upstream edges of features, -100 means start 100 bp upstream of first base of feature, 0 means start exactly at first base of feature
				  #for downstream edges of feature, -99 means start 99 bp upstream of last base of feature, 1 means start exactly after last base of feature
bin_end = 1999  #Distance from end position of last bin to edge of feature. "-" means upstream. 
			  #e.g., for upstream edges of features, -1 means end exactly before first base of feature. 99 means end 99 bp downstream of first base of feature.
		      #for downstream edges of feature 0 means end at last base of feature, 100 means end 100 bp downstream of last base of feature.
				
how_many_chr = 12

import glob
import os
import statistics



#Set up to run through every input file in directory
for filename in glob.glob(partial_ATCGmap_path + '/' + '*' + ending):

	#I/O
	ATCGmap_path = filename
	out_path = filename[:-15] + output_ending
	out_file = open(out_path, 'w') #Write-to file

	#####################################################################################################################################
	#####################################################################################################################################
	#Define a function Bins_Meth_Analyzer that calculates average methylation for each set of bins that are derived from the input GFF file
	######################################################################################################################################
	#####################################################################################################################################
	def Bins_Meth_Analyzer(i):
	
		################################################################################################################
		#INITIALIZE SOME MORE VARIABLES
		################################################################################################################
	
		#I/O
		ATCGmap_path = filename
		ATCGmap_file = open(ATCGmap_path) #methylated genome file
		ATCGmap_file_is_open = True
		GFF_file = open(GFF_path) #File with coordinates in GFF format
		GFF_file_is_open = True
		
		#list of rightmost feature end for each chromosome. This is important for identifying overlapping features that will need to be stored as "residuals"
		ends = [0 for x in range(how_many_chr)]  

		#list of residual features
		residuals = [] 
		residual_count = 0
		residual_number = 0

		#Initialize an empty list to hold the methyl counts for each bp of a feature
		#This list will hold the following:
		#     0 - percent mC
		#	  1 - context--CG, CHG, CHH, CCH, or other
		#	  2 - coverage, as measured by number of reads that overlap with this position
		single_bp_values = [0,'',0]

		total_bin_length = 0 #for calculating coverage later
		
		#Initialize some variables to hold summary values
		CG_list = []
		CHG_list = []
		CHH_list = []
		CCH_list = []
		covered_bp_counts = 0
		read_nt_counts = 0
 
		#Initialize an empty list to hold the final summary values for all features 
		#This list will hold the following: 
		#     0 - relative percent methylated CG
		#	  1 - relative percent methylate CHG 
		#	  2 - relative percent methylated CHH
		#	  3 - relative percent methylated CCH	
		#	  4 - number of CG's covered by at least one read in the feature 
		#	  5 - number of CHG's covered by at least one read in the feature
		#     6 - number of CHH's covered by at least one read in the feature
		#	  7 - number of CCH's covered by at least one read in the feature	
		#     8 - coverage, as measured by percent of bp in region that are covered by at least one read (0-1)
		#	  9 - coverage, as measured by number of nt in reads divided by length of region ( >= 0)
		whole_feature_values = [0 for value in range(10)] 

		#########################################################################################################
		#Define a function, "GFF_line_extractor" for getting information from a GFF line
		##########################################################################################################
		
		def GFF_line_extractor(residual_count, residual_number, total_bin_length):
			
			#For the first run through the map file/GFF file
			if GFF_file_is_open:		
			
				#get details on current feature (but skip it if overlaps with a prior one) 
				skip = True #assume this feature overlaps one that was already processed and should be skipped for now and added to residuals list
				while skip:	
				
					GFF_line = GFF_file.readline() #Get the next feature

					#Get info from GFF line
					if GFF_line != '': #file.readline() yields '' at end of file
												
						GFF_cols = GFF_line.split('\t')
						GFF_chr = int(GFF_cols[0]) - 1
				
						if upstream_edge:
							if GFF_cols[6] == '+':
								GFF_start = int(GFF_cols[3]) + (bin_start + i*bin_size)
							else:
								GFF_start = int(GFF_cols[4]) - (bin_start + (i+1)*bin_size - 1)
						else:
							if GFF_cols[6] == '+':
								GFF_start  = int(GFF_cols[4]) + (bin_start + i*bin_size)
							else:
								GFF_start = int(GFF_cols[3]) - (bin_start + (i+1)*bin_size - 1)
	
						GFF_end = GFF_start + bin_size - 1		
						total_bin_length += bin_size # keep track of total length for coverage calculation later

						#Keep track of feature ends and add feature to residuals if there is an overlap.
						#Then get a new feature
						if GFF_start <= ends[GFF_chr] and GFF_start > 0:
							residuals.append(GFF_line) #add overlapping feature to end of residuals list
							residual_count += 1
		
						else:
							ends[GFF_chr] = GFF_end
							skip = False
							map_ahead_of_GFF = True #don't get a new map position
					
					#if the end of the GFF file is reached, need to skip through the rest of the map file
					else:	
						GFF_chr = how_many_chr #This number is 1 higher than map_chr - 1 can be. It will cause all subsequent map positions to be skipped
						map_ahead_of_GFF = False #time to get a new map position
						return GFF_chr, '', '', residual_count, residual_number, map_ahead_of_GFF, total_bin_length
						print("you should not be reading this")
				
				
				#in special case where ATCGmap ends before GFF_file, need to run through the rest of the GFF_file to record that additional features exist and that they have no coverage
				#this is not necessary later when processing residuals because these features will not be identified as residuals, even if they overlap. 
				if ATCGmap_file_is_open == False:
					for GFF_line in GFF_file:
						total_bin_length += bin_size # keep track of total length for coverage calculation later	
		
			#For the subsequent runs through the map file/residuals list
			else:
				
				#get details on current feature (but skip it if overlaps with a prior one) 
				skip = True #assume this feature overlaps one that was already processed and should be skipped for now and added to residuals list
				while skip:
				
					if residual_number < old_residual_count: #need to watch out for end of current residuals list
						GFF_line = residuals[residual_number] #get the next feature
						residual_number += 1 #keep track of which one to get next time
					
						GFF_cols = GFF_line.split('\t')
						GFF_chr = int(GFF_cols[0]) - 1
				
						if upstream_edge:
							if GFF_cols[6] == '+':
								GFF_start = int(GFF_cols[3]) + (bin_start + i*bin_size)
							else:
								GFF_start = int(GFF_cols[4]) - (bin_start + (i+1)*bin_size - 1)
						else:
							if GFF_cols[6] == '+':
								GFF_start  = int(GFF_cols[4]) + (bin_start + i*bin_size)
							else:
								GFF_start = int(GFF_cols[3]) - (bin_start + (i+1)*bin_size - 1)
	
						GFF_end = GFF_start + bin_size - 1		

						#Keep track of feature ends and add feature to residuals if there is an overlap.
						#Then get a new feature
						if GFF_start <= ends[GFF_chr] and GFF_start > 0:
							residuals.append(GFF_line) #add overlapping feature to end of residuals list
							residual_count += 1
		
						else:
							ends[GFF_chr] = GFF_end
							skip = False
							map_ahead_of_GFF = True #don't get a new map position
							
					
					#if the last current residual is reached, need to skip through the rest of the map file
					else:
						GFF_chr = how_many_chr #This number is 1 higher than map_chr - 1 can be. It will cause all subsequent map positions to be skipped
						map_ahead_of_GFF = False #time to get a new map position
						return GFF_chr, '', '', residual_count, residual_number, map_ahead_of_GFF, total_bin_length
						print("you should not be reading this")
						
			return GFF_chr, GFF_start, GFF_end, residual_count, residual_number, map_ahead_of_GFF, total_bin_length
			print("you should not be reading this")

		#########################################################################################################
		#Run through ATCGmap file in parallel with GFF file, analyzing features one by one
		##########################################################################################################

		#get info from first GFF line
		GFF_chr, GFF_start, GFF_end, residual_count, residual_number, map_ahead_of_GFF, total_bin_length = GFF_line_extractor(residual_count, residual_number, total_bin_length)

		#get next map position
		for map_line in ATCGmap_file:	
			map_cols = map_line.split('\t')
			try:
				map_chr = int(map_cols[0]) - 1 
				numeric_chr = True
			except ValueError:
				numeric_chr = False
				
			if numeric_chr:	
				map_position = int(map_cols[2])	

				#Whenever map position is ahead of feature end position, call in a new feature from GFF file or residuals list using GFF_line_extractor function.  
				map_ahead_of_GFF = True		
				while map_ahead_of_GFF: #In this while loop, new features will be called until 
											#A) the feature is determined to be ahead of the map position
											#B) the map position is determined to be within the feature
											#C) there are no more features in the GFF file			
		
					#Option 1: map chromosome is ahead 
					if GFF_chr < map_chr: 
						# get another GFF feature
						GFF_chr, GFF_start, GFF_end, residual_count, residual_number, map_ahead_of_GFF, total_bin_length = GFF_line_extractor(residual_count, residual_number, total_bin_length)
						
					#Option 2: GFF chromosome is ahead
					elif GFF_chr > map_chr:
						map_ahead_of_GFF = False #break out of while loop and get next map position
			 
					#Option 3: GFF chromosome matches map chromosome
					else:

						#Option 3.1: map position to the right of feature
						if map_position > GFF_end: 
				
							# get another GFF feature
							GFF_chr, GFF_start, GFF_end, residual_count, residual_number, map_ahead_of_GFF, total_bin_length = GFF_line_extractor(residual_count, residual_number, total_bin_length)

						#Option 3.2: map position to the left of feature
						elif map_position < GFF_start:
							map_ahead_of_GFF = False #break out of while loop and get next map position
							
						#Option 3.3: map position is in GFF feature 
						else:
						
							#add methylation information for current map position to single_bp_values
							
							covered_bp_counts += 1 #Count the number of bp covered by a read
							read_nt_counts += int(map_cols[5]) + int(map_cols[6]) + int(map_cols[7]) + int(map_cols[8]) + int(map_cols[10]) + int(map_cols[11]) + int(map_cols[12]) + int(map_cols[13]) 
							
							try:
								single_bp_values[0] = float(map_cols[15])
								
								if map_cols[3] == 'CHH':
									CHH_list.append(single_bp_values[0]) #Record number of CHHs in region and the methylation frequency of each 
									if map_cols[4] == 'CC':
										CCH_list.append(single_bp_values[0]) #Record number of CCHs in region and the methylation frequency of each 
								elif map_cols[3] == 'CG':
									CG_list.append(single_bp_values[0]) #Record number of CGs in region and the methylation frequency of each 
								elif map_cols[3] == 'CHG':
									CHG_list.append(single_bp_values[0]) #Record number of CHGs in region and the methylation frequency of each 										
		
							except ValueError:
								pass
							map_ahead_of_GFF = False #break out of while loop and get next map position (it doesn't matter whether map position is actually ahead)
	
	
		
		ATCGmap_file.close()
		ATCGmap_file_is_open = False		
	
		#Run through rest of GFF file if necessary
		GFF_chr, GFF_start, GFF_end, residual_count, residual_number, map_ahead_of_GFF, total_bin_length = GFF_line_extractor(residual_count, residual_number, total_bin_length)
		
		GFF_file.close()
		GFF_file_is_open = False
		#########################################################################################################
		#Run through ATCGmap file and residual features in parallel, processing features one by one
		##########################################################################################################
	 

		#Run through ATCGmap file to analyze residuals (overlapping features that were skipped)
		#Repeat this until there are no more residuals. 
		while residual_count > 0:
			print("now analyzing " + str(residual_count) + " overlapping features that were skipped in the last pass")
			del(residuals)[:residual_number] #Delete already-analyzed residuals
	
			old_residual_count = residual_count #to prevent recycling residuals on the current pass through map file
			residual_number = 0 #tracks position within residual list on current pass through map file 
			residual_count = 0 #counts residuals added to list on current pass through map file
	
			#reopen map file
			ATCGmap_file = open(ATCGmap_path) #methylated genome file
			ATCGmap_file_is_open = True
		
			#list of rightmost feature end for each chromosome. This is important for identifying overlapping features that will need to be stored as residuals
			ends = [0 for x in range(how_many_chr)]  
	
			#get info from first GFF line
			GFF_chr, GFF_start, GFF_end, residual_count, residual_number, map_ahead_of_GFF, total_bin_length = GFF_line_extractor(residual_count, residual_number, total_bin_length)
	
			#get next map position
			for map_line in ATCGmap_file:	
				map_cols = map_line.split('\t')
				try:
					map_chr = int(map_cols[0]) - 1 
					numeric_chr = True
				except ValueError:
					numeric_chr = False
				
				if numeric_chr:	
					map_position = int(map_cols[2])	

					#Whenever map position is ahead of feature end position, call in a new feature from GFF file or residuals list using GFF_line_extractor function.  
					map_ahead_of_GFF = True		
					while map_ahead_of_GFF: #In this while loop, new features will be called until 
												#A) the feature is determined to be ahead of the map position
												#B) the map position is determined to be within the feature
												#C) there are no more features in the GFF file			

														
						#Option 1: map chromosome is ahead 
						if GFF_chr < map_chr: 				
							GFF_chr, GFF_start, GFF_end, residual_count, residual_number, map_ahead_of_GFF, total_bin_length = GFF_line_extractor(residual_count, residual_number, total_bin_length)
				
						#Option 2: GFF chromosome is ahead
						elif GFF_chr > map_chr:
							map_ahead_of_GFF = False #break out of while loop and get next map position
		 
						#Option 3: matching chromosomes
						else:
	
							#Option 3.1: map position is to the right of feature
							if map_position > GFF_end: 
								GFF_chr, GFF_start, GFF_end, residual_count, residual_number, map_ahead_of_GFF, total_bin_length = GFF_line_extractor(residual_count, residual_number, total_bin_length)
						
							#Option 3.2: map position is to the left of feature
							elif map_position < GFF_start:
								map_ahead_of_GFF = False #break out of while loop and get next map position
				
							#Option 3.3: map position is in GFF feature 
							else:
								#add methylation information for current map position to single_bp_values
								
								covered_bp_counts += 1 #Count the number of bp covered by a read
								read_nt_counts += int(map_cols[5]) + int(map_cols[6]) + int(map_cols[7]) + int(map_cols[8]) + int(map_cols[10]) + int(map_cols[11]) + int(map_cols[12]) + int(map_cols[13]) 
								
								try:
									single_bp_values[0] = float(map_cols[15])
									
									if map_cols[3] == 'CHH':
										CHH_list.append(single_bp_values[0]) #Record number of CHHs in region and the methylation frequency of each 
										if map_cols[4] == 'CC':
											CCH_list.append(single_bp_values[0]) #Record number of CCHs in region and the methylation frequency of each
									elif map_cols[3] == 'CG':
										CG_list.append(single_bp_values[0]) #Record number of CGs in region and the methylation frequency of each 
									elif map_cols[3] == 'CHG':
										CHG_list.append(single_bp_values[0]) #Record number of CHGs in region and the methylation frequency of each 										
									
								except ValueError:
									pass
								map_ahead_of_GFF = False #break out of while loop and get next map position (it doesn't matter whether map position is actually ahead)							
					
			#close map file			
			ATCGmap_file.close()
	
	
		##################################################################################################
		#Calculate means, etc
		##################################################################################################
	
		#Get relative fraction methylated CG
		try:
			whole_feature_values[0] = statistics.mean(CG_list)	
		except statistics.StatisticsError:
			whole_feature_values[0] = ''

		#Get relative fraction methylated CHG
		try:
			whole_feature_values[1] = statistics.mean(CHG_list)	
		except statistics.StatisticsError:
			whole_feature_values[1] = ''
	
		#Get relative fraction methylated CHH
		try:
			whole_feature_values[2] = statistics.mean(CHH_list)	
		except statistics.StatisticsError:
			whole_feature_values[2] = ''
			
		#Get relative fraction methylated CCH
		try:
			whole_feature_values[3] = statistics.mean(CCH_list)	
		except statistics.StatisticsError:
			whole_feature_values[3] = ''	
 
		'''#Record number of CGs in the feature that are covered by at least one read
		whole_feature_values[3] = len(CG_list)

		#Record number of CHGs in the feature that are covered by at least one read
		whole_feature_values[4] = len(CHG_list)

		#Record number of CHHs in the feature that are covered by at least one read
		whole_feature_values[5] = len(CHH_list)'''
	
		print("total length of bin " + str(i + 1) + ": " + str(total_bin_length))
		#Calculate and record the fraction of the feature that is covered by at least one read
		whole_feature_values[8] = covered_bp_counts/total_bin_length

		#Calculate and record the coverage (total nt in reads over length of region) 
		whole_feature_values[9] = read_nt_counts/total_bin_length

		'''#Calculate and record absolute percent methylated CG
		if whole_feature_values[0] != '':
			whole_feature_values[8] = whole_feature_values[0]*CG_count/covered_nt_counts
		else:
			whole_feature_values[8] = ''
	
		#Calculate and record absolute percent methylated CHG
		if whole_feature_values[1] != '':
			whole_feature_values[9] = whole_feature_values[1]*CHG_count/covered_nt_counts
		else:
			whole_feature_values[9] = ''
	
		#Calculate and record absolute percent methylated CHH
		if whole_feature_values[2] != '':
			whole_feature_values[10] = whole_feature_values[2]*CHH_count/covered_nt_counts
		else:
			whole_feature_values[10] = '':
		'''

		#WRITE OUTPUT (THIS COULD BE MODIFIED USING ABOVE DETAILS AND STATISTICS MODULE TO INCLUDE STANDARD ERRORS FOR EACH)
		out_file.write('\t' + str(whole_feature_values[0]) + '\t' + str(whole_feature_values[1]) + '\t' + str(whole_feature_values[2]) + '\t' + str(whole_feature_values[3]) + '\t' + str(whole_feature_values[8]) + '\t' + str(whole_feature_values[9]) + '\n')


	#################################################################################################################
	#################################################################################################################
	#Call on Bin_Meth_Analyzer function for each bin
	################################################################################################################
	#################################################################################################################

	out_file.write('bin start\tCG\tCHG\tCHH\tCCH\tcoverage (fraction of length)\tcoverage\n')
	for i in range((bin_end + 1 - bin_start)//bin_size):
		out_file.write(str(bin_start + i*bin_size))
		print("working on bin " + str(i + 1) + " of " + str(len(range((bin_end + 1 - bin_start)//bin_size))))
		
		Bins_Meth_Analyzer(i)


	################################################################################################################
	#Close I/O
	################################################################################################################

	out_file.close()
