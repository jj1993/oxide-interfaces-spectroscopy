'''
Created on 13 April 2014

@author: dqz93389 Christoph
'''
from utils.dRangeUtil import drange
from scan.concurrentAnalyserScan import pathscan
from functions.functionClassFor2Scannables import ScannableFunctionClassFor2Scannables
from functions.ratioscannable import ratioFun

psn_list=[]

def realign_Jbranch():
    fsj1('In')
    old_exitslit_y=ss4ygap.getPosition()
    old_exitslit_x=ss4xgap.getPosition()
    ss4xgap.asynchronousMoveTo(0.2) 
    ss4ygap.asynchronousMoveTo(0.2)
    while ss4xgap.isBusy() or ss4ygap.isBusy():
        sleep(0.1)
        
    rscan sm3fpitch -2 2 0.1 sm5iamp8
    pos sm3fpitch peak.result.pos
    
    ss4xgap.asynchronousMoveTo(old_exitslit_x) 
    ss4ygap.asynchronousMoveTo(old_exitslit_y)
    while ss4xgap.isBusy() or ss4ygap.isBusy():
        sleep(0.1)
    print()
    print('Sm3fpitch realigned. Exit slit at old position!')
    
    
    
def i_switch_hard():
    E_1storder = ienergy()
    E_3rdorder = E_1storder*3
    igap.moveTo(ienergy.calc(E_3rdorder, 5))
    rscan igap -0.3 0.3 0.02 hm3iamp20
    pos igap peak.result.pos
    print('Switched to', E_3rdorder, 'keV.')
    print('Remember to measure in kinetic energies!')

def i_switch_soft():
    E_1storder = ienergy()
    igap.moveTo(ienergy.calc(ienergy(),3))
    

def save_psn(psn_r):
	psn_list=[]
	if psn_r-1 > len(psn_list):
		print psn_r, 'is to high! Use ', len(psn_list),' instead to create a new entry.'
	elif psn_r-1 < len(psn_list) and psn_r > 0:
		print 'Overwriting ', psn_r,'.'
		psn_list[psn_r-1]=[smpmx.getPosition(),smpmy.getPosition(),smpmz.getPosition(),smpmpolar.getPosition()]
		print 'smpmx: ', psn_list[psn_r-1][0]
		print 'smpmy: ', psn_list[psn_r-1][1]
		print 'smpmz: ', psn_list[psn_r-1][2]
		print 'smpmpolar: ', psn_list[psn_r-1][3]
	elif psn_r-1 == len(psn_list):
		print 'Position saved as position', psn_r,'.'
		psn_list.append([smpmx.getPosition(),smpmy.getPosition(),smpmz.getPosition(),smpmpolar.getPosition()])
		print 'smpmx: ', psn_list[psn_r-1][0]
		print 'smpmy: ', psn_list[psn_r-1][1]
		print 'smpmz: ', psn_list[psn_r-1][2]
		print 'smpmpolar: ', psn_list[psn_r-1][3]

def movetopsn(psn_r):
	if psn_r-1 >= len(psn_list):
		print 'Position', pos_nr, 'is not yet defined!'
	elif psn_r-1 < len(psn_list):
		smpmx.asynchronousMoveTo(psn_list[psn_r-1][0])
		smpmy.asynchronousMoveTo(psn_list[psn_r-1][1])
		smpmz.asynchronousMoveTo(psn_list[psn_r-1][2])
		smpmpolar.asynchronousMoveTo(psn_list[psn_r-1][3])
		while smpmx.isBusy() or smpmy.isBusy() or smpmz.isBusy() or smpmpolar.isBusy():
			sleep(0.1)
		val=psn_r-1
		print 'position' , psn_r, ' reached:'
		print 'smpmx: ' , psn_list[psn_r-1][0]
		print 'smpmy: ' , psn_list[psn_r-1][1]
		print 'smpmz: ' , psn_list[psn_r-1][2]
		print 'smpmpolar: ' , psn_list[psn_r-1][3]
		
#del ratio39to8
#ratio39to8 = ScannableFunctionClassFor2Scannables("ratio39to8", "smpmiamp39", "sm5iamp8", ratioFun)


def NEXAFS_variable(Energy1,Energy2,EnergyStep,preedge,preedge_step,postedge,postedge_step):
	#def NEXAFS_variable(Energy1,Energy2,EnergyStep,preedge,preedge_step,postedge,postedge_step):
	from utils.dRangeUtil import drange
	from scan.concurrentAnalyserScan import pathscan
	
	# Select scan type:
	# scantype 1= ResPES
	# scantype 2= NEXAFS
	# scantype 3= NEXAFS with undulator offset
	scan_type=3
	gap_offset=0.0
	
	
	#####################  
	# This constructs a set of photon energies
	#####################  
	myenergyPath= drange(Energy1-preedge,Energy1-preedge_step,preedge_step) + drange(Energy1,Energy2,EnergyStep) + drange(Energy2+postedge_step,Energy2+postedge,postedge_step)
	jgapPath=[]
	for n in range(0,len(myenergyPath)):
		jgapPath.append(jenergy.calc(myenergyPath[n],1)-gap_offset)
		myenergyPath[n]=myenergyPath[n]*1000
		
	myenergyPath_zip = zip(jgapPath,myenergyPath)
	myenergyPath_list = []
	for l in myenergyPath_zip:
	     myenergyPath_list.append(list(l))
	     myenergyPath_tuple = tuple(myenergyPath_list)
	
	##################### 
	# Select and start scan
	#####################     
	
	fsj1('Out')
	oldtime=time.time()
	
	if scan_type==1:
	    analyserscan((jgap,pgmenergy), myenergyPath_tuple, ew4000, "user.seq", sm5amp8,1, smpmamp39,1) 
	if scan_type==2:
	    analyserscan((jenergy,y), myenergyPath_tuple, sm5amp8,.5, smpmamp39,.5) 
	if scan_type==3:
	    #analyserscan((jgap,pgmenergy), myenergyPath_tuple, sm5iamp8, smpmiamp39) 
	    #analyserscan((jgap,pgmenergy), myenergyPath_tuple, sm5iamp8, smpmiamp39) #smpmamp39
	    analyserscan((jgap,pgmenergy), myenergyPath_tuple, sm5amp8, 1, smpmamp39, 1, smpm)
	newtime=time.time()
	print 'Scan took:', newtime-oldtime, 'seconds'
	#fsj1('In')
	
	
	##############   Test with dummy motors
	#analyserscan((x, y), myenergyPath_tuple, nixswr, hm3iamp20, smpmamp39)   # this will measure nixswReflectivity only
	#analyserscan((x, y), myenergyPath_tuple, ew4000, "user.seq",nixswr, hm3iamp20, smpmamp39)  # this will measure XPS and nixswReflectivity

        
#def HAXPES_scans():
#    analyserscan ew4000 "Survey.seq"
#    print 'Survey scan finished'

    
    
#def SXPS_scans():
#    analyserscan ew4000 "Survey_NCA_SX.seq"
#    print 'Survey scan finished'




def VO2_XSW():
	y_start=-552
	Ebragg=4.1974

	pos smpmy y_start
	for i in range(0,1):
	    interruptable()
	    print '--- starting iteration number: ', i+1, 'out of 1 ---'
	    pos smpmy y_start+i*1
	    pos fsi1 'Out'
	    scan dcmenergy Ebragg-0.0015 Ebragg+0.0015 0.0001 hm3iamp20 lakeshoreC3 nixswr
	    Ebragg = peak.result.pos
	    analyserscan dcmenergy Ebragg-0.0015 Ebragg+0.0015 0.0001 ew4000 "XSW_fixed.seq" hm3amp20 2 smpm lakeshoreC3 nixswr
	    pos fsi1 'In'   
	print 'XSW scans finished' 
    
def Julias_scans():
    energies= [ 0.875,
               0.835,
               0.72,
               0.65,
               0.475,
               0.325,
               0.273,
               0.245]
    
    sequences = ['Corelevels_soft_LTO_LMO.seq',
                 'Corelevels_soft_Mn2p.seq',
                 'Corelevels_soft_O1s.seq',
                 'Corelevels_soft_Ti2p.seq',
                 'Corelevels_soft_C1s.seq',
                 'Corelevels_soft_P2p.seq',
                 'Corelevels_soft_Mn3s.seq',
                 'Corelevels_soft_Li1s.seq']
    
    i=0
    smpmy1 = smpmy.getPosition()
    smpmy_step = 0.05
    
    oldtime=time.time()
    
    for energy in energies:
        interruptable()
        print i, energy, sequences[i]
        pos jenergy energy
        analyserscan ew4000 sequences[i] sm5iamp8 smpm; pos fsj1 'In'
        pos smpmy smpmy1-smpmy_step*i
        i=i+1
    
    newtime=time.time()
    print 'Scan took:', (newtime-oldtime)/60, 'minutes'
    
    pos smpmy smpmy1-smpmy_step*len(energies)
    pos ienergy 2.35
    pos igap 14
    analyserscan ew4000 'Corelevels_2350_LTO_LMO.seq' hm3iamp20 smpm; pos fsi1 'In'
    
    pos smpmy smpmy1-smpmy_step*(len(energies)+1)
    pos ienergy 2.35
    pos igap 8.83
    analyserscan ew4000 'Corelevels_7050_LTO_LMO.seq' hm3iamp20 smpm; pos fsi1 'In'
    
    print 'Measurements finished!'


def ienergy_off(target):
    pos fsi1 'In'
    gap_offset=0.1
    dcmenergy.moveTo(target)
    igap_off=ienergy.calc(target,3)+gap_offset
    igap.moveTo(igap_off)

def rockingcurves(total,center):
    for i in range(0,total):
        interruptable()
        scan dcmenergy center-0.0015 center+0.0015 0.0001 nixswr hm3amp20 1 smpmamp39 1 nixswr
        print '--- iterration number: ', i+1, 'out of', total, '---'
        
        
def battery_scans(smpmy1):
    # start scans 7keV
    smpmy_step=0.05
    nstep = 4
    pos igap 8.83
    for i in range(0,nstep):
            interruptable()
            print '--- iteration number: ', i+1
            analyserscan ew4000 "Corelevels_7050_b.seq" smpm
            pos smpmy smpmy1-smpmy_step*i
    analyserscan ew4000 "Survey_7050_b.seq" smpm
    pos fsi1 'In'
    print 'Sample 7keV finished'
     
    pos fsj1 'In'