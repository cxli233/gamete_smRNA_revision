#This uses the information in a BS Seeker 2 ATCGmap file to determine the mean methylation level in every feature in a GFF file.
#It also produces the mean, median, standard deviation and standard error for the whole set of features. 
#GFF file must be sorted numerically, chromosome than start position. 
#Overlapping GFFs are tolerated as long as they are sorted.
#GFF file should only include numeric chromosomes, no scaffolds or organellar.
#ATCGmap file must be sorted numerically, chromosome than position (sort -n -k1,1 -k3,3). 
#All non-numeric chromosomes are skipped. 




#OUTLINE
#For each feature, determine methylation status at each bp and record in single_bp_values.
	#Use GFF_line_extractor function to extract details form each feature
#When a feature is completed, calculate average values for the feature using whole_feature_summarizer function. 
#Then save the summary values in all_feature_summaries and write them to out_file.
#When all features have been analyzed, calculate average of all features and print the values. 
#If there are no overlapping features, the above will occur in a single pass through ATCGmap file. 
#If there are overlapping features, multiple runs through will be made the ATCGmap file, each time with a smaller set of features 
     #(the "residuals" or overlapping feature that were missed in the prior run)
#Repeat above for every input file in directory

################################################################################################################
#INITIALIZE SOME VARIABLES
################################################################################################################


#I/O

#Directory of ATCGmaps (unlike previous versions, don't include '/' at end!)
partial_ATCGmap_path = '/scratch/gent/rice_BS/sorted_maps'

#ATCGmap file endings
ending = ".ATCGmap.sorted"

#GFF file
GFF_path = '/home/gent/references/sperm-specific_24siRNAs.cov' #File with coordinates in GFF format

#output ending
output_ending =  ending.split('.')[0] + GFF_path.split('/')[-1][:GFF_path.split('/')[-1].rfind('.')] + '.meth'

how_many_chr = 12

import glob
import os
import statistics

#summary file
summaryFile = open(partial_ATCGmap_path + '/' + output_ending[:-5] + '.BGFsummary.txt', 'w')
summaryFile.write('\tmCG\tmCHG\tmCHH\tmCCH\n')

#Set up to run through every input file in directory
for filename in glob.glob(partial_ATCGmap_path + '/' + '*' + ending):

	################################################################################################################
	#INITIALIZE SOME MORE VARIABLES
	################################################################################################################
	
	#I/O
	ATCGmap_path = filename
	out_path = filename[:-len(ending)] + '_' + output_ending
	ATCGmap_file = open(ATCGmap_path) #methylated genome file
	GFF_file = open(GFF_path) #File with coordinates in GFF format
	GFF_file_is_open = True
	out_file = open(out_path, 'w') #Write-to file

	#list of rightmost feature end for each chromosome. This is important for identifying overlapping features that will need to be stored as "residuals"
	ends = [0 for x in range(how_many_chr)]  

	#list of residual features
	residuals = [] 
	residual_count = 0
	residual_number = 0
	covered_bp_counts = 0
	read_nt_counts = 0
	
	#Initialize an empty list to hold the final summary values for each feature 
	#This list will hold the following: 
	#     0 - relative percent methylated CG
	#	  1 - relative percent methylate CHG 
	#	  2 - relative percent methylated CHH
	#	  3 - relative percent methylated CHH
	#	  4 - number of CG's covered by at least one read in the feature 
	#	  5 - number of CHG's covered by at least one read in the feature
	#     6 - number of CHH's covered by at least one read in the feature
	#	  7 - number of CCH's covered by at least one read in the feature	
	#     8 - coverage, as measured by percent of bp in region that are covered by at least one read (0-1)
	#	  9 - coverage, as measured by number of nt in reads divided by length of region ( >= 0)
	whole_feature_values = [0 for value in range(10)] 

	#Initialize an empty list to hold the final summary values for all features
	#This list will hold the following:
	#     0 - relative percent methylated CG
	#	  1 - relative percent methylated CHG 
	#	  2 - relative percent methylated CHH
	#	  3 - relative percent methylated CCH
	#     4 - read coverage	
	all_features_summary = [[],[],[],[],[]]

	out_file.write('chr\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tmCG\tmCHG\tmCHH\tmCCH\tcovered CGs\tcovered CHGs\tcovered CHHs\tcovered CCHs\tcovered bases\tcoverage\n')

	#########################################################################################################
	#Define a function, "GFF_line_extractor" for getting information from a GFF line and resetting single_bp_values
	##########################################################################################################

	def GFF_line_extractor(covered_bp_counts, read_nt_counts, first_feature_already_processed, residual_count, residual_number):
		
		#invoke whole_feature_summarizer on previous feature
		if first_feature_already_processed:
			whole_feature_summarizer(covered_bp_counts, read_nt_counts)
		
		#For the first run through the map file/GFF file
		if GFF_file_is_open:		
		
			#get details on current feature (but skip it if overlaps with a prior one) 
			skip = True
			while skip:	
				
				GFF_line = GFF_file.readline() #Get the next feature	
				
				#Get info from GFF line
				if GFF_line != '': #file.readline() yields '' at end of file
											
					GFF_cols = GFF_line.split('\t')
					GFF_chr = int(GFF_cols[0]) - 1
					GFF_start = int(GFF_cols[3])
					GFF_end = int(GFF_cols[4])
	
					#Keep track of feature ends and add feature to residuals if there is an overlap.
					#Then get a new feature
					if GFF_start <= ends[GFF_chr] and GFF_start > 0:
						residuals.append(GFF_line) #add overlapping feature to end of residuals list
						residual_count += 1
	
					else:
						ends[GFF_chr] = GFF_end
						skip = False
						map_ahead_of_GFF = True #don't get a new map position
						out_file.write(GFF_line[:-1]) #output all info from GFF file still in GFF format
						
				#if the end of the GFF file is reached, need to skip through the rest of the map file
				else:	
					GFF_chr = how_many_chr #This number is 1 higher than map_chr - 1 can be. It will cause all subsequent map positions to be skipped
					map_ahead_of_GFF = False #time to get a new map position
					return GFF_chr, '', '', [], residual_count, residual_number, map_ahead_of_GFF, 0, 0 
					print("you should not be reading this")
				
			#in special case where ATCGmap ends before GFF_file, need to run through the rest of the GFF_file to record that additional features exist and that they have no coverage
			#this is not necessary later when processing residuals because these features will not be identified as residuals, even if they overlap. 
			if ATCGmap_file_is_open == False:
				for GFF_line in GFF_file:
					out_file.write(GFF_line[:-1])
					whole_feature_summarizer(covered_bp_counts, read_nt_counts)
				
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
					GFF_start = int(GFF_cols[3])
					GFF_end = int(GFF_cols[4])
		
					#Keep track of feature ends and add feature to residuals if there is an overlap.
					#Then get a new feature
					if GFF_start <= ends[GFF_chr] and GFF_start > 0:
						residuals.append(GFF_line) #add overlapping feature to end of residuals listfread_nt
						residual_count += 1
	
					else:
						ends[GFF_chr] = GFF_end
						skip = False
						map_ahead_of_GFF = True #don't get a new map position
						out_file.write(GFF_line[:-1]) #output all info from GFF file still in GFF format
				
				#if the last current residual is reached, need to skip through the rest of the map file
				else:
					GFF_chr = how_many_chr #This number is 1 higher than map_chr - 1 can be. It will cause all subsequent map positions to be skipped
					map_ahead_of_GFF = False #time to get a new map position
					return GFF_chr, '', '', [], residual_count, residual_number, map_ahead_of_GFF, 0, 0
					print("you should not be reading this")
						
		#Initialize an empty list to hold the methyl counts for each bp of a feature
		#This list will hold the following:
		#     0 - percent mCG
		#	  1 - context--CG, CHG, CHH, CCH, or other (xyz- format)
		#	  2 - coverage, as measured by number of reads that overlap with this position
		single_bp_values = [[0,'',0] for bp in range(GFF_end - GFF_start + 1)]
		covered_bp_counts = 0
		read_nt_counts = 0
			
		return GFF_chr, GFF_start, GFF_end, single_bp_values, residual_count, residual_number, map_ahead_of_GFF, covered_bp_counts, read_nt_counts
		print("you should not be reading this")		

		

	#########################################################################################################
	#Define a function, "whole_feature_summarizer" for calculating methylation averages for a whole feature
	##########################################################################################################
	#The parameter single_bp_values contains a list of methyl counts corresponding to a single feature

	def whole_feature_summarizer(covered_bp_counts, read_nt_counts):
		
		if len(single_bp_values) > 0:
			#Initialize lists to hold methylation frequencies for each covered CG, CHG, and CHH position and to hold various counts
			CG_list = []
			CHG_list = []
			CHH_list = []	
			CCH_list = []		
	

			#Go through feature, one nt at a time, and count methyls
			for bp in range(len(single_bp_values)):
			
				if single_bp_values[bp][1] == 'CHH':
					CHH_list.append(single_bp_values[bp][0]) #Record number of CHHs in region and the methylation frequency of each
				elif single_bp_values[bp][1] == 'CCH':
					CCH_list.append(single_bp_values[bp][0]) #Record number of CCHs in region and the methylation frequency of each 
					CHH_list.append(single_bp_values[bp][0]) #Record number of CHHs in region and the methylation frequency of each
				elif single_bp_values[bp][1] == 'CG':
					CG_list.append(single_bp_values[bp][0]) #Record number of CGs in region and the methylation frequency of each 
				elif single_bp_values[bp][1] == 'CHG':
					CHG_list.append(single_bp_values[bp][0]) #Record number of CHGs in region and the methylation frequency of each 	
	
		
			#Calculate summary values
	
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
	 
			#Record number of CGs in the feature that are covered by at least one read
			whole_feature_values[4] = len(CG_list)
	
			#Record number of CHGs in the feature that are covered by at least one read
			whole_feature_values[5] = len(CHG_list)
	
			#Record number of CHHs in the feature that are covered by at least one read
			whole_feature_values[6] = len(CHH_list)
			
			#Record number of CCHs in the feature that are covered by at least one read
			whole_feature_values[7] = len(CCH_list)

			#Calculate and record the fraction of the feature that is covered by at least one read (not separated by strand, so maximum value is 1)
			whole_feature_values[8] = covered_bp_counts/(GFF_end - GFF_start + 1)

			#Calculate and record the coverage (total nt in reads over length of region) 
			whole_feature_values[9] = read_nt_counts/(GFF_end - GFF_start + 1)
	
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
				whole_feature_values[10] = '' '''
		
			#write summary to the list "all_features_summary" and to output file
			for i in range(4):
				out_file.write('\t' + str(whole_feature_values[i])) #write mCG, mCHG, mCHH, mCCH frequencies to out_file
				if whole_feature_values[i] != '':
					all_features_summary[i].append(whole_feature_values[i])	#add mCG, mCHG, mCHH, mCCH frequencies to all_features_summary
			all_features_summary[4].append(read_nt_counts) #add coverage to all_features_summary	
			for i in range(4,10):	
				out_file.write('\t' + str(whole_feature_values[i])) #write other details to out_file
			out_file.write('\n')
		
	#########################################################################################################
	#Run through GFF file and ATCGmap file in parallel, analyzing features one by one
	##########################################################################################################
	first_feature_already_processed = False
	ATCGmap_file_is_open = True
	
	#Get GFF_feature chromosome and position and reset single_bp_values

	GFF_chr, GFF_start, GFF_end, single_bp_values, residual_count, residual_number, map_ahead_of_GFF, covered_bp_counts, read_nt_counts = GFF_line_extractor(covered_bp_counts, read_nt_counts, first_feature_already_processed, residual_count, residual_number)
	first_feature_already_processed = True #this is necessary to prevent invokation of whole_feature_summarizer before the first feature has been processed

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
					GFF_chr, GFF_start, GFF_end, single_bp_values, residual_count, residual_number, map_ahead_of_GFF, covered_bp_counts, read_nt_counts = GFF_line_extractor(covered_bp_counts, read_nt_counts, first_feature_already_processed, residual_count, residual_number)
					
				#Option 2: GFF chromosome is ahead
				elif GFF_chr > map_chr:
					map_ahead_of_GFF = False #get next map position
			 
				#Option 3: GFF chromosome matches map chromosome
				else:
		
					#Option 3.1: map position to the right of feature
					if map_position > GFF_end: 
									
						# get another GFF feature
						GFF_chr, GFF_start, GFF_end, single_bp_values, residual_count, residual_number, map_ahead_of_GFF, covered_bp_counts, read_nt_counts = GFF_line_extractor(covered_bp_counts, read_nt_counts, first_feature_already_processed, residual_count, residual_number)
					
					#Option 3.2: map position to the left of feature
					elif map_position < GFF_start:
						map_ahead_of_GFF = False #get next map position
					
					#Option 3.3: map position is in GFF feature 
					else:
						#add methylation information for current map position to single_bp_values
						bp = map_position - GFF_start
						covered_bp_counts += 1   #Count the number of bp covered by a read					
						single_bp_values[bp][1] = map_cols[3]
						if map_cols[3] == 'CHH':
							if map_cols[4] == 'CC':
								single_bp_values[bp][1] = 'CCH'
						read_nt_counts += int(map_cols[5]) + int(map_cols[6]) + int(map_cols[7]) + int(map_cols[8]) + int(map_cols[10]) + int(map_cols[11]) + int(map_cols[12]) + int(map_cols[13]) 

						try:
							single_bp_values[bp][0] = float(map_cols[15])
						except ValueError:
							pass
						map_ahead_of_GFF = False #get next map position
	
	
	ATCGmap_file.close()
	ATCGmap_file_is_open = False		
	
	#Call whole_feature_summarizer function on any remaining features in GFF file
	GFF_chr, GFF_start, GFF_end, single_bp_values, residual_count, residual_number, map_ahead_of_GFF, covered_bp_counts, read_nt_counts = GFF_line_extractor(covered_bp_counts, read_nt_counts, first_feature_already_processed, residual_count, residual_number)
	
	GFF_file.close()
	GFF_file_is_open = False

	#########################################################################################################
	#Run through residuals and ATCGmap file in parallel, processing features one by one
	##########################################################################################################
	 

	#Run through ATCGmap file to analyze residuals (overlapping features that were skipped)
	#Repeat this until there are no more residuals. 
	while residual_count > 0:
		print("now analyzing " + str(residual_count) + " overlapping features that were skipped in the last pass")
		del(residuals)[:residual_number] #Delete already-analyzed residuals
	
		old_residual_count = residual_count #to prevent recycling residuals on the current pass through map file
		residual_number = 0 #tracks position within residual list on current pass through map file 
		residual_count = 0 #counts residuals added to list on current pas through map file
	
		first_feature_already_processed = False
	
		#reopen map file
		ATCGmap_file = open(ATCGmap_path) #methylated genome file
		ATCGmap_file_is_open = True
		
		#list of rightmost feature end for each chromosome. This is important for identifying overlapping features that will need to be stored as residuals
		ends = [0 for x in range(how_many_chr)]  
	
		#Get first GFF_feature chromosome and position and reset single_bp_values
		GFF_chr, GFF_start, GFF_end, single_bp_values, residual_count, residual_number, map_ahead_of_GFF, covered_bp_counts, read_nt_counts = GFF_line_extractor(covered_bp_counts, read_nt_counts, first_feature_already_processed, residual_count, residual_number)
		first_feature_already_processed = True #this is necessary to prevent invokation of whole_feature_summarizer before the first feature has been processed
	
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
						GFF_chr, GFF_start, GFF_end, single_bp_values, residual_count, residual_number, map_ahead_of_GFF, covered_bp_counts, read_nt_counts = GFF_line_extractor(covered_bp_counts, read_nt_counts, first_feature_already_processed, residual_count, residual_number)
					
					#Option 2: GFF chromosome is ahead
					elif GFF_chr > map_chr:
						map_ahead_of_GFF = False #get next map position
			 
					#Option 3: matching chromosomes
					else:
		
						#Option 3.1: map position is to the right of feature
						if map_position > GFF_end: 
							GFF_chr, GFF_start, GFF_end, single_bp_values, residual_count, residual_number, map_ahead_of_GFF, covered_bp_counts, read_nt_counts = GFF_line_extractor(covered_bp_counts, read_nt_counts, first_feature_already_processed, residual_count, residual_number)

						#Option 3.2: map position is to the left of feature
						elif map_position < GFF_start:
							map_ahead_of_GFF = False #get next map position
					
						#Option 3.3: map position is in GFF feature 
						else:
							#add methylation information for current map position to single_bp_values
							bp = map_position - GFF_start
							covered_bp_counts += 1   #Count the number of bp covered by a read		
							try:
								single_bp_values[bp][0] = float(map_cols[15])
								single_bp_values[bp][1] = map_cols[3]
								if map_cols[3] == 'CHH':
									if map_cols[4] == 'CC':
										single_bp_values[bp][1] = 'CCH'
								read_nt_counts += int(map_cols[5]) + int(map_cols[6]) + int(map_cols[7]) + int(map_cols[8]) + int(map_cols[10]) + int(map_cols[11]) + int(map_cols[12]) + int(map_cols[13]) 
							except ValueError:
								pass
							map_ahead_of_GFF = False #get next map position
		
		ATCGmap_file.close()
		ATCGmap_file_is_open = False		
	
		#Call whole_feature_summarizer for last feature
		GFF_chr, GFF_start, GFF_end, single_bp_values, residual_count, residual_number, map_ahead_of_GFF, covered_bp_counts, read_nt_counts = GFF_line_extractor(covered_bp_counts, read_nt_counts, first_feature_already_processed, residual_count, residual_number)
	
		GFF_file.close()
		GFF_file_is_open = False

	##########################################################################################################
	#PRINT SUMMARY
	##########################################################################################################
	
	summaryFile.write(ATCGmap_path[ATCGmap_path.rfind('/') + 1:] + '\t' + str(statistics.mean(all_features_summary[0])) + '\t' + str(statistics.mean(all_features_summary[1])) + '\t' + str(statistics.mean(all_features_summary[2])) + '\t' + str(statistics.mean(all_features_summary[3])) + '\n')
	print("File: " + out_path)
	print("CG mean: " + str(statistics.mean(all_features_summary[0])))
	print("CHG mean: " + str(statistics.mean(all_features_summary[1])))
	print("CHH mean: " + str(statistics.mean(all_features_summary[2])))
	print("CCH mean: " + str(statistics.mean(all_features_summary[3])))
	print("CG standard deviation: " + str(statistics.stdev(all_features_summary[0])))
	print("CHG standard deviation: " + str(statistics.stdev(all_features_summary[1])))
	print("CHH standard deviation: " + str(statistics.stdev(all_features_summary[2])))
	print("CCH standard deviation: " + str(statistics.stdev(all_features_summary[3])))
	print("CG standard error of the mean: " + str(statistics.stdev(all_features_summary[0])/len(all_features_summary[0])))
	print("CHG standard error of the mean: " + str(statistics.stdev(all_features_summary[1])/len(all_features_summary[1])))
	print("CHH standard error of the mean: " + str(statistics.stdev(all_features_summary[2])/len(all_features_summary[2])))
	print("CCH standard error of the mean: " + str(statistics.stdev(all_features_summary[3])/len(all_features_summary[3])))
	print()

	out_file.close()
summaryFile.close()