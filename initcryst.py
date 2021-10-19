import peakopt
import runcrystfel
import crystplots

#Peak search optimization
#peakopt.main()

#Unit cell parameters distribution, indexing, integration, merging, calculation of figures of merit using indexamajig and conversions to mtz or xscale files
runcrystfel.main()

#Comparatevely plot of figures of merit in resolution shells for different subdatasets (labels) in different paths (dir)
#label=[['16'],['32'],['48'],['64']]
#dir=['../20210611_16/compare_hkl_part_liso_11_0/','../20210616_32/compare_hkl_part_liso_11_0/','../20210625_48/compare_hkl_part_liso_11_0/', '../20210627_64/compare_hkl_part_liso_11_0/']
#dir=['../20210611_16/check_hkl_part_liso_11_/','../20210616_32/check_hkl_part_liso_11_/','../20210625_48/check_hkl_part_liso_11_/', '../20210627_64/check_hkl_part_liso_11_/']

#crystplots.plot_check(label,'part','0', dir)
#crystplots.plot_compare(label,'part','0', dir)
