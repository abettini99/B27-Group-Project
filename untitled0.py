# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 13:16:52 2020

@author: natas
"""

#<<<<<<< HEAD
#motion = 'spiral'    # set motion - 'phugoid', 'shortperiod', 'aperroll', 'dutchroll', 'dutchrollYD', 'spiral'
#data   = manouvre(data, motion)                          # sliced data array for phugoid motion
#=======
#
#<<<<<<< HEAD
#if motion in ['aperroll', 'dutchroll', 'dutchrollYD', 'spiral']:
#    sysa = ctl.StateSpace(Aa, Ba, C, D)                                                  # create state-space system for symmetric eigenmotions
#    evals, evacs = eig(As)                                                               # compute eigenvalues and eigenvectors
#    print('=================== EIGENVALUES OF MATRIX Aa ==================')
#    print(evals)
#    print('===============================================================')
#
#    tstop = data.time.iloc[-1] - data.time.iloc[0]                                       # normalise final time value for manouvre
#    dt  = np.arange(0, tstop + 0.1, 0.1)                                                 # create time vector with 0.1s step size
#
#    units = ['[rad]', '[rad]', '[rad/s]', '[rad/s]']                                             # list with units of columns for plotting
#    u = [np.radians(data.delta_a), np.radians(data.delta_r)]                             # [rad] input array given input at each time for [da, dr]
#    columns = [r'\beta', r'\phi', r'p', r'r']                    # names of invidiual columns for DataFrame
#    eigenmotion = []                                                                     # initialise empty list
#
#    flightdata = [np.radians(data.Ahrs1_Roll), np.radians(data.Ahrs1_bRollRate), np.radians(data.Ahrs1_bYawRate)]
#
#    t, y, x = ctl.forced_response(sysa, dt, U=u)                                         # calculate forced response
#    df2 = pd.DataFrame(np.transpose(y), columns=columns)                                 # convert forced response to DataFrame
#    eigenmotion.append(df2)                                                              # append DataFrame to individual list
#    eigenmotion = pd.concat(eigenmotion, axis=1)                                         # concatenate list into panda dataframe along axis 1
#
#    if motion == 'aperroll':
#        fig1, ax1 = plt.subplots(4,1, squeeze=False, figsize=(16,9))                     # initialise figure with 4 rows and 1 column
#        for i in range(0,3):
#            ax1[i,0].plot(t, eigenmotion.iloc[:,i+1], 'C1', label='Numerical Model')     # plot each variable from output vector
#            ax1[i,0].plot(t, flightdata[i], c='k', label='Experimental Data')            # plot each variable from test flight data
#            ax1[i,0].set_xlabel('$t$ [s]')                                               # set label of x-axis
#            ax1[i,0].set_ylabel('${}$ {}'.format(eigenmotion.columns[i+1], units[i+1]))  # set label of y-axis
#            ax1[i,0].minorticks_on()                                                     # set minor ticks
#            ax1[i,0].grid(which='major', linestyle='-', linewidth='0.5', color='black')  # customise major grid
#            ax1[i,0].grid(which='minor', linestyle=':', linewidth='0.5', color='grey')   # customise minor grid
#            ax1[i,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')        # set legend for subplot
#
#        ax1[3,0].plot(t, u[0], c='k', linestyle='--', label='Aileron Deflection')        # plot input variable delta_a
#        ax1[3,0].plot(t, u[1], c='k', linestyle='-',label='Rudder Deflection')           # plot input variable delta_r
#        ax1[3,0].set_xlabel('$t$ [s]')                                                   # set label of y-axis
#        ax1[3,0].set_ylabel('$\delta_a, \delta_r$ [rad]')                                # set label of y-axis
#        ax1[3,0].minorticks_on()                                                         # set minor ticks
#        ax1[3,0].grid(which='major', linestyle='-', linewidth='0.5', color='black')      # customise major grid
#        ax1[3,0].grid(which='minor', linestyle=':', linewidth='0.5', color='grey')       # customise minor grid
#        ax1[3,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')            # set legend
#
#        fig1.tight_layout(pad=1.0)                                                       # increase spacing between subplots
#        fig1.suptitle('Aperiodic Roll')                                                  # set title of figure
#        # fig1.savefig('images/aperiodicroll.png', dpi=300, bbox_inches='tight')           # save figure
#
#    if motion == 'dutchroll':
#        fig1, ax1 = plt.subplots(4,1, squeeze=False, figsize=(16,9))                     # initialise figure with 4 rows and 1 column
#        for i in range(0,3):
#            ax1[i,0].plot(t, eigenmotion.iloc[:,i+1], 'C1', label='Numerical Model')     # plot each variable from output vector
#            ax1[i,0].plot(t, flightdata[i], c='k', label='Experimental Data')            # plot each variable from test flight data
#            ax1[i,0].set_xlabel('$t$ [s]')                                               # set label of x-axis
#            ax1[i,0].set_ylabel('${}$ {}'.format(eigenmotion.columns[i+1], units[i+1]))  # set label of y-axis
#            ax1[i,0].minorticks_on()                                                     # set minor ticks
#            ax1[i,0].grid(which='major', linestyle='-', linewidth='0.5', color='black')  # customise major grid
#            ax1[i,0].grid(which='minor', linestyle=':', linewidth='0.5', color='grey')   # customise minor grid
#            ax1[i,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')        # set legend for subplot
#
#        ax1[3,0].plot(t, u[0], c='k', linestyle='--', label='Aileron Deflection')        # plot input variable delta_a
#        ax1[3,0].plot(t, u[1], c='k', linestyle='-',label='Rudder Deflection')           # plot input variable delta_r
#        ax1[3,0].set_xlabel('$t$ [s]')                                                   # set label of y-axis
#        ax1[3,0].set_ylabel('$\delta_a, \delta_r$ [rad]')                                # set label of y-axis
#        ax1[3,0].minorticks_on()                                                         # set minor ticks
#        ax1[3,0].grid(which='major', linestyle='-', linewidth='0.5', color='black')      # customise major grid
#        ax1[3,0].grid(which='minor', linestyle=':', linewidth='0.5', color='grey')       # customise minor grid
#        ax1[3,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')            # set legend
#
#        fig1.tight_layout(pad=1.0)                                                       # increase spacing between subplots
#        fig1.suptitle('Dutch Roll without Yaw Damper')                                   # set title of figure
#        # fig1.savefig('images/dutchroll.png', dpi=300, bbox_inches='tight')               # save figure
#
#    if motion == 'dutchrollYD':
#        fig1, ax1 = plt.subplots(4,1, squeeze=False, figsize=(16,9))                     # initialise figure with 4 rows and 1 column
#        for i in range(0,3):
#            ax1[i,0].plot(t, eigenmotion.iloc[:,i+1], 'C1', label='Numerical Model')     # plot each variable from output vector
#            ax1[i,0].plot(t, flightdata[i], c='k', label='Experimental Data')            # plot each variable from test flight data
#            ax1[i,0].set_xlabel('$t$ [s]')                                               # set label of x-axis
#            ax1[i,0].set_ylabel('${}$ {}'.format(eigenmotion.columns[i+1], units[i+1]))  # set label of y-axis
#            ax1[i,0].minorticks_on()                                                     # set minor ticks
#            ax1[i,0].grid(which='major', linestyle='-', linewidth='0.5', color='black')  # customise major grid
#            ax1[i,0].grid(which='minor', linestyle=':', linewidth='0.5', color='grey')   # customise minor grid
#            ax1[i,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')        # set legend for subplot
#
#        ax1[3,0].plot(t, u[0], c='k', linestyle='--', label='Aileron Deflection')        # plot input variable delta_a
#        ax1[3,0].plot(t, u[1], c='k', linestyle='-',label='Rudder Deflection')           # plot input variable delta_r
#        ax1[3,0].set_xlabel('$t$ [s]')                                                   # set label of y-axis
#        ax1[3,0].set_ylabel('$\delta_a, \delta_r$ [rad]')                                # set label of y-axis
#        ax1[3,0].minorticks_on()                                                         # set minor ticks
#        ax1[3,0].grid(which='major', linestyle='-', linewidth='0.5', color='black')      # customise major grid
#        ax1[3,0].grid(which='minor', linestyle=':', linewidth='0.5', color='grey')       # customise minor grid
#        ax1[3,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')            # set legend
#
#        fig1.tight_layout(pad=1.0)                                                       # increase spacing between subplots
#        fig1.suptitle('Dutch Roll with Yaw Damper')                                      # set title of figure
#        # fig1.savefig('images/dutchrollYD.png', dpi=300, bbox_inches='tight')             # save figure
#
#    if motion == 'spiral':
#        fig1, ax1 = plt.subplots(4,1, squeeze=False, figsize=(16,9))                     # initialise figure with 4 rows and 1 column
#        for i in range(0,3):
#            ax1[i,0].plot(t, eigenmotion.iloc[:,i+1], 'C1', label='Numerical Model')     # plot each variable from output vector
#            ax1[i,0].plot(t, flightdata[i], c='k', label='Experimental Data')            # plot each variable from test flight data
#            ax1[i,0].set_xlabel('$t$ [s]')                                               # set label of x-axis
#            ax1[i,0].set_ylabel('${}$ {}'.format(eigenmotion.columns[i+1], units[i+1]))  # set label of y-axis
#            ax1[i,0].minorticks_on()                                                     # set minor ticks
#            ax1[i,0].grid(which='major', linestyle='-', linewidth='0.5', color='black')  # customise major grid
#            ax1[i,0].grid(which='minor', linestyle=':', linewidth='0.5', color='grey')   # customise minor grid
#            ax1[i,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')        # set legend for subplot
#
#        ax1[3,0].plot(t, u[0], c='k', linestyle='--', label='Aileron Deflection')        # plot input variable delta_a
#        ax1[3,0].plot(t, u[1], c='k', linestyle='-',label='Rudder Deflection')           # plot input variable delta_r
#        ax1[3,0].set_xlabel('$t$ [s]')                                                   # set label of y-axis
#        ax1[3,0].set_ylabel('$\delta_a, \delta_r$ [rad]')                                # set label of y-axis
#        ax1[3,0].minorticks_on()                                                         # set minor ticks
#        ax1[3,0].grid(which='major', linestyle='-', linewidth='0.5', color='black')      # customise major grid
#        ax1[3,0].grid(which='minor', linestyle=':', linewidth='0.5', color='grey')       # customise minor grid
#        ax1[3,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')            # set legend
#
#        fig1.tight_layout(pad=1.0)                                                       # increase spacing between subplots
#        fig1.suptitle('Spiral')                                                          # set title of figure
#        # fig1.savefig('images/spiralroll.png', dpi=300, bbox_inches='tight')              # save figure
#
#plt.show()
##
##data = importdata('flightdata.mat')
##new_data = manouvre(data, 'aperroll')
##new_data.plot(x='time', y='Ahrs1_bRollRate')
#=======