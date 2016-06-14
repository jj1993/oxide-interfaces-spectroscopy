"""
FUNCTIONS TO DO GOLD CALIBRATION MEASUREMENTS FOR THE UvA GROUP
"""

def do_gold_scanning():
    """
    Script that does measurements for the UvA group on a plain gold sample
    for calibration purposes.
    """

##    # Position the sample
##    pos smpmx (...)
##    pos smpmy (...)
##    pos smpmz (...)
##    pos smpmpolar (...)
##    pos smpmazimuth (...)
##    while (smpmx.isBusy() or smpmy.isBusy() or smpmz.isBusy()
##           or smpmpolar.isBusy() or smpmazimuth.isBusy()):
##        sleep(0.1)

    realign_Jbranch()
    realign_Ibranch()

    while True:
        # Continues for infinity
        try: megamindy()
        except: print "ALL SOFT SWEEP SERIES FAILED.....???"
        try: megadeep()
        except: print "ALL SOFT DEPT SCANS FAILED.....???"
        try: hard()
        except: print "ALL HARD MEASUREMENTS FAILED.....???"
        try: soft()
        except: print "ALL SOFT MEASUREMENTS FAILED.....???"

    return

def megamindy():
    ################################
    ################################
    # SrMnSb_2 Megamindy sweep map #
    ################################
    ################################

    our_ss4xgap = 0.03
    hv = .36

    jenergy.setPolarisation('LH')
    pos ss4xgap our_ss4xgap
    pos jenergy hv
    while ss4xgap.isBusy() or jenergy.isBusy():
        sleep(0.1)

    #### NOTE
    #### SrMnSb.seq must have:
    # Swept
    # pass energy 70
    # Kinetic energies 349-356 (:/20)
    # size 20meV

    print '====================================='
    print
    print 'SrMnSb_2 MEGAMINDY SWEEP SERIES'
    print 'Measurement nr:\t',n+1
    print
    print '====================================='

    sequence = 'megasweep.seq'

    try:
        # Do the acual scan!
        analyserscan ew4000 sequence
        pos fsj1 'In'
    except:
        print "SCAN FAILED. CONTINUES ON NEXT SCAN!"
        continue

    return

def megadeep():
    ###########################
    ###########################
    # SrMnSb_2 Megamindy dept #
    ###########################
    ###########################

    our_ss4xgap = 0.03
    sequence = 'megadepth.seq'

    jenergy.setPolarisation('LH')
    pos ss4xgap our_ss4xgap
    while ss4xgap.isBusy():
        sleep(0.1)

    #### NOTE
    #### SrMnSb.seq must have:
    # Swept
    # pass energy 70
    # Kinetic energies 349-356 (:/20)
    # size 20meV
    
    for n, t in enumerate([.22,.24,.26,.28,.30,.32,.34,.36,.38,.40,.42,.44,.46,.48,.50]):
        interruptable()
        
        print '====================================='
        print
        print 'SrMnSb_2 MEGAMINDY DEPT SCANS'
        print 'Measurement nr:\t',n+1
        print
        print '====================================='

        try:
            hv = t
            
            pos jenergy hv
            while jenergy.isBusy():
                sleep(0.1)
                
            # Do the acual scan!
            analyserscan ew4000 sequence
            pos fsj1 'In'
        except:
            print "SCAN FAILED. CONTINUES ON NEXT SCAN!"
            continue

    return


def hard():
    hard_measurements = [  #Energy     igap         pass-energy
                            #(...,       ...,           ...)
                            (2.35,       14.261,      'hardlowpass100.seq'),
                            (2.35,       14.261,      'hardlowpass70.seq'),
                            (3.0,       22.1264,        'hardlowpass100.seq'),
                            (3.0,       22.1264,        'hardlowpass70.seq')
                            #(2.35,      8.839 ,     'hardhighpass200.seq') # THIS IS THE 7.05keV SETTING!!
                        ]

    for n, m in enumerate(hard_measurements):
        interruptable()
        
        # New variables are defined
        energy = m[0]
        our_igap = m[1]
        sequence = m[2]
        
        print '====================================='
        print
        print 'HARD MEASUREMENT
        print 'Measurement nr:\t',n+1
        print 'Energy:\t\t\t',energy
        print 'Slid width:\t\t',our_igap
        print 'Sequence:\t\t',sequence
        print
        print '====================================='
            
        # Select new slidth and energy settings
        pos ienergy energy
        while ienergy.isBusy():
            sleep(0.1)
        pos igap our_igap
        while igap.isBusy():
            sleep(0.1)            

        # Do the acual scan! Each sequence does two sweeps
        analyserscan ew4000 sequence
        pos fsj1 'In'

    return

def soft():

    # FIRST 3 ENERGIES
    for energy in [.45915, .46015, .46285]:

        pos jenergy energy
        while jenergy.isBusy():
            sleep(0.1)

        print '====================================='
        print
        print 'SOFT MEASUREMENT
        print 'Energy:\t\t\t',energy
        print
        print '====================================='

        # 3 sweeps + 3 fixed for 3 energies, total of 18 files
        for t in [('fermiswept70.seq',0.06),('fermifixed70.seq',0.06),
                  ('fermiswept100.seq',0.06),('fermifixed100.seq',0.06)]:
            sequence, our_ss4x = t
            
            ##### DETERMINE ss4x GAP AT 3 PASS ENERGIES!! ######
            pos ss4xgap our_ss4xgap
            while ss4xgap.isBusy():
                sleep(0.1)
            ##### DETERMINE ss4x GAP AT 3 PASS ENERGIES!! ######
                
            # Do the acual scan!
            analyserscan ew4000 sequence
            pos fsj1 'In'
            
    # 4TH ENERGY  
    energy = .642
    pos jenergy energy
    while jenergy.isBusy():
        sleep(0.1)

    print '====================================='
    print
    print 'SOFT MEASUREMENT
    print 'Energy:\t\t\t',energy
    print
    print '====================================='

    # total of 6 sweeps
    for t in [('au3pswept20.seq',0.02),
              ('fermiswept70.seq',0.06),('au3pswept70.seq',0.06),
              ('fermiswept100.seq',0.06),('au3pswept100.seq',0.06)]:
        sequence, our_ss4x = t
        
        ##### DETERMINE ss4x GAP AT 3 PASS ENERGIES!! ######
        pos ss4xgap our_ss4xgap
        while ss4xgap.isBusy():
            sleep(0.1)
        ##### DETERMINE ss4x GAP AT 3 PASS ENERGIES!! ######
            
        # Do the acual scan!
        analyserscan ew4000 sequence
        pos fsj1 'In'

    # 5TH ENERGY
    energy = .900      
    pos jenergy energy
    while jenergy.isBusy():
        sleep(0.1)

    print '====================================='
    print
    print 'SOFT MEASUREMENT
    print 'Energy:\t\t\t',energy
    print
    print '====================================='

    # total of 6 sweeps
    for t in [('au3pswept20.seq',0.02),
              ('fermiswept70.seq',0.06),('au3pswept70.seq',0.06),
              ('fermiswept100.seq',0.06),('au3pswept100.seq',0.06)]:
        sequence, our_ss4x = t
        
        ##### DETERMINE ss4x GAP AT 3 PASS ENERGIES!! ######
        pos ss4xgap our_ss4xgap
        while ss4xgap.isBusy():
            sleep(0.1)
        ##### DETERMINE ss4x GAP AT 3 PASS ENERGIES!! ######
            
        # Do the acual scan!
        analyserscan ew4000 sequence
        pos fsj1 'In'

    return

"""
END OF UVA FUNCTIONS
"""
