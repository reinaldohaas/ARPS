#################################################################################
# FILE:		plttimeseries.py
# CREATED:      2011-06-29 by Brett Roberts
# MODIFIED:     2011-07-01 by Brett Roberts
# REQUIRES:	matplotlib
#
# Produces timeseries plots for various quantitites along trajectories.
# The quantities themselves are specified within the script, and can be
# modified to suit the user's needs. For each quantitity (or set) specified,
# a separate plot will be created for each trajectory in the input file.
#
# *** TO CUSTOMIZE WHICH VARIABLES ARE PLOTTED, PLEASE EDIT THE SECTION
#     OF THE CODE INDICATED BELOW ***
#
#
# Usage:
#       python /abspath/plttimeseries.py runname.trajc_xxx output_dir
#
#       NOTES:
#       - You must explicitly call python at the command line, rather than
#          simply executing the script (e.g., ./plttimeseries.py will not work)
#       - You must specify the absolute path of this script at the command
#          line, even in the event that you are already in the directory
#          which contains it
################################################################################

import sys
import os
import matplotlib.pyplot as plt
import matplotlib as mpl

filepath = os.path.dirname(sys.argv[0])
sys.path.append(filepath + '/../modules/')

# Import arpscalctrajc library
import arpscalctrajc_lib as actl

# Number of command-line arguments
n_args = len(sys.argv)

# Check arguments
if(n_args != 3):
        print "ERROR: You must specify an arpscalctrajc file as input, as well as an output directory."
        exit
else:
    calcfile = os.path.abspath(sys.argv[1])
    output_dir = os.path.abspath(sys.argv[2])

# ************* CODE FOR USER TO EDIT BEGINS HERE ***************

#################################################################
# TIMESERIES TO BE PLOTTED					#
#################################################################

# Base series to extract from the calctrajc file
# [ label, series ]
# *** You must include quantities you will sum, integrate, or differentiate ***
base_series = [
 [ 'z' ,		'z' ],
 [ 'w' ,		'w' ],
 [ 'u' ,                'u' ],
 [ 'v' ,                'v' ],
 [ 'ptprt',		'ptprt' ], 
 [ 'vortx' ,		'vortx' ],
 [ 'vorty' ,		'vorty' ],
 [ 'vortz' ,		'vortz' ],
 [ 'vorts' ,            'vorts' ],
 [ 'vortc' ,            'vortc' ],
 [ 'vortx_{gen}',	'vortx_gen'],
 [ 'vortx_{stch}',	'vortx_stch'],
 [ 'vortx_{tilt}',	'vortx_tilt'],
 [ 'vortx_{mix}',      	'vortx_mix'],
 [ 'vorty_{gen}',	'vorty_gen'],
 [ 'vorty_{stch}',      'vorty_stch'],
 [ 'vorty_{tilt}',      'vorty_tilt'],
 [ 'vorty_{mix}',       'vorty_mix'],
 [ 'vortz_{stch}',      'vortz_stch'],
 [ 'vortz_{tilt}',      'vortz_tilt'],
 [ 'vortz_{mix}',	'vortz_mix'],
 [ 'vorts_{gen}',       'vorts_gen'],
 [ 'vorts_{stch}',      'vorts_stch'],
 [ 'vorts_{tilt}',      'vorts_tilt'],
 [ 'vorts_{mix}',       'vorts_mix'],
 [ 'vorts_{exhg}',      'vorts_exhg'],
 [ 'vortc_{gen}',       'vortc_gen'],
 [ 'vortc_{stch}',      'vortc_stch'],
 [ 'vortc_{tilt}',      'vortc_tilt'],
 [ 'vortc_{mix}',       'vortc_mix'],
 [ 'vortc_{exhg}',      'vortc_exhg'],
]

# Series comprised of the sum of any set of series above
# [ label, [series1, series2, ...] ]
summed_series = [
 [ 'vortx_{prog-dt}' ,		['vortx_{gen}','vortx_{stch}','vortx_{tilt}','vortx_{mix}'] ],
 [ 'vorty_{prog-dt}' ,    	['vorty_{gen}','vorty_{stch}','vorty_{tilt}','vorty_{mix}'] ],
 [ 'vortx_{prog-dt-nomix}' ,	['vortx_{gen}','vortx_{stch}','vortx_{tilt}'] ],
 [ 'vorty_{prog-dt-nomix}' , 	['vorty_{gen}','vorty_{stch}','vorty_{tilt}'] ],
 [ 'vortz_{prog-dt}' ,    	['vortz_{stch}','vortz_{tilt}','vortz_{mix}'] ],
 [ 'vorts_{prog-dt}' ,          ['vorts_{gen}','vorts_{stch}','vorts_{tilt}','vorts_{exhg}','vorts_{mix}'] ],
 [ 'vortc_{prog-dt}' ,          ['vortc_{gen}','vortc_{stch}','vortc_{tilt}','vortc_{exhg}','vortc_{mix}'] ],
 [ 'vorts_{prog-dt-nomix}' ,    ['vorts_{gen}','vorts_{stch}','vorts_{tilt}','vorts_{exhg}'] ],
 [ 'vortc_{prog-dt-nomix}' ,    ['vortc_{gen}','vortc_{stch}','vortc_{tilt}','vortc_{exhg}'] ],
]

# Time derivative of any series above
# [ label, series ]
differentiated_series = [
 [ 'vortx_{dt}',	'vortx' ],
 [ 'vorty_{dt}',        'vorty' ],
 [ 'vortz_{dt}',	'vortz' ],
 [ 'vorts_{dt}',        'vorts' ],
 [ 'vortc_{dt}',        'vortc' ],
]

# Time integral of any series above
# [ label, series, init_series]
integrated_series = [
 [ 'vortx_{prog}',		'vortx_{prog-dt}',		'vortx' ],
 [ 'vorty_{prog}',		'vorty_{prog-dt}',		'vorty' ],
 [ 'vortx_{prog-nomix}',      	'vortx_{prog-dt-nomix}',      	'vortx' ],
 [ 'vorty_{prog-nomix}',      	'vorty_{prog-dt-nomix}',      	'vorty' ],
 [ 'vortz_{prog}',		'vortz_{prog-dt}',		'vortz' ],
 [ 'vorts_{prog}',		'vorts_{prog-dt}',		'vorts' ],
 [ 'vortc_{prog}',		'vortc_{prog-dt}',      	'vortc' ],
 [ 'vorts_{prog-nomix}',	'vorts_{prog-dt-nomix}',	'vorts' ],
 [ 'vortc_{prog-nomix}',	'vortc_{prog-dt-nomix}',	'vortc' ],
]

# Difference between any series above
# [label, [base_series, subtract_series] ]
difference_series = [
 [ 'vortx_{dt-err}',	['vortx_{prog-dt}','vortx_{dt}'] ],
 [ 'vorty_{dt-err}',    ['vorty_{prog-dt}','vorty_{dt}'] ],
 [ 'vortz_{dt-err}',    ['vortz_{prog-dt}','vortz_{dt}'] ],
 [ 'vorts_{dt-err}',    ['vorts_{prog-dt}','vorts_{dt}'] ],
 [ 'vortc_{dt-err}',    ['vortc_{prog-dt}','vortc_{dt}'] ],
]

# Relative error using any two of the series above
# [label, [base_series, comp_series] ]
relerror_series = [
 [ 'vortx_{dt-relerr}',    ['vortx_{dt}','vortx_{prog-dt}'] ],
 [ 'vorty_{dt-relerr}',    ['vorty_{dt}','vorty_{prog-dt}'] ],
]

#################################################################
# PLOTS TO CREATE						#
#################################################################

# [ [label1, label2, ...], units, yaxisstr, titlestr, filename ]
plots = [
 [ ['ptprt'],                                   'K',            'ptprt',                'ptprt',                'ptprt' ],
 [ ['u'],                                       'm \ \ s^{-1}', 'u',                    'u',                    'u' ],
 [ ['v'],                                       'm \ \ s^{-1}', 'v',                    'v',                    'v' ],
 [ ['w'],                                       'm \ \ s^{-1}', 'w',                    'w',                    'w' ],
 [ ['z'],                                       'm',            'z',                    'z',                    'z' ],
 [ ['vortx_{prog}','vortx_{prog-nomix}','vortx'],                       's^{-1}',       '\omega_{x}',           '\omega_{x}',           'vortx' ],
 [ ['vorty_{prog}','vorty_{prog-nomix}','vorty'],                       's^{-1}',       '\omega_{y}',           '\omega_{y}',           'vorty' ],
 [ ['vorts_{prog}','vorts_{prog-nomix}','vorts'],                       's^{-1}',       '\omega_{s}',           '\omega_{s}',           'vorts' ],
 [ ['vortc_{prog}','vortc_{prog-nomix}','vortc'],                       's^{-1}',       '\omega_{c}',           '\omega_{c}',           'vortc' ],
 [ ['vortz_{prog}','vortz'],                    's^{-1}',       '\omega_{z}',           '\omega_{z}',           'vortz' ],
 [ ['vortx_{dt-err}'],                          's^{-2}',       '\omega_{x_{err}}',     '\omega_{x_{err}}',     'vortxerr' ],
 [ ['vorty_{dt-err}'],                          's^{-2}',       '\omega_{y_{err}}',     '\omega_{y_{err}}',     'vortyerr' ],
 [ ['vorts_{dt-err}'],                          's^{-2}',       '\omega_{s_{err}}',     '\omega_{s_{err}}',     'vortserr' ],
 [ ['vortc_{dt-err}'],                          's^{-2}',       '\omega_{c_{err}}',     '\omega_{c_{err}}',     'vortcerr' ],
 [ ['vortx_{dt-relerr}'],                          's^{-2}',       '\omega_{x_{relerr}}',  '\omega_{x_{relerr}}',  'vortxrelerr' ],
 [ ['vorty_{dt-relerr}'],                          's^{-2}',       '\omega_{y_{relerr}}',  '\omega_{y_{relerr}}',  'vortyrelerr' ],
 [ ['vortx_{dt}','vortx_{prog-dt}'],            's^{-2}',       'd{\omega}_{x}/dt',     'd{\omega}_{x}/dt',     'vortx_dt' ],
 [ ['vorty_{dt}','vorty_{prog-dt}'],            's^{-2}',       'd{\omega}_{y}/dt',     'd{\omega}_{y}/dt',     'vorty_dt' ],
 [ ['vorts_{dt}','vorts_{prog-dt}'],            's^{-2}',       'd{\omega}_{s}/dt',     'd{\omega}_{s}/dt',     'vorts_dt' ],
 [ ['vortc_{dt}','vortc_{prog-dt}'],            's^{-2}',       'd{\omega}_{c}/dt',     'd{\omega}_{c}/dt',     'vortc_dt' ],
 [ ['vortx_{gen}','vortx_{stch}','vortx_{tilt}','vortx_{mix}'], 's^{-2}',       'd\omega_{x}/dt',       'd\omega_{x}/dt',       'vortxcpts' ],
 [ ['vorty_{gen}','vorty_{stch}','vorty_{tilt}','vorty_{mix}'], 's^{-2}',       'd\omega_{y}/dt',       'd\omega_{y}/dt',       'vortycpts' ],
 [ ['vorts_{gen}','vorts_{stch}','vorts_{tilt}','vorts_{mix}','vorts_{exhg}'], 's^{-2}',       'd\omega_{s}/dt',       'd\omega_{s}/dt',       'vortscpts' ],
 [ ['vortc_{gen}','vortc_{stch}','vortc_{tilt}','vortc_{mix}','vortc_{exhg}'], 's^{-2}',       'd\omega_{c}/dt',       'd\omega_{c}/dt',       'vortccpts' ],
]

# ************** CODE FOR USER TO EDIT ENDS HERE ****************

# Open calctrajc file
try:
	calcdata = actl.get_trajc_array(calcfile)
except IOError:
        print "ERROR: The input file you specified either is not a valid arpstrajc file or does not exist."
        exit

header = actl.get_header(calcfile)
ntimes = actl.get_npnts(calcfile)
ntrajcs = actl.get_ntrajcs(calcfile)
nvars = len(header)
run_name = actl.get_run_name(calcfile)

mylib = actl.calctrajcdata(calcfile)

# Manually set drawing parameters
mpl.rcParams['axes.linewidth'] = 0.5
mpl.rcParams['mathtext.default'] = 'sf'
legendfont = mpl.font_manager.FontProperties(size=9)

# Produce a single timeseries plot
def calcplot( filename, times, timeseries_array, label_array, title, xlabel, ylabel ):
	print "\nCreating plot " + title + "..."

	fig = plt.figure(num=None, figsize=(10,6), dpi=96, facecolor='w')

	ax = fig.add_subplot(111)

	ax.grid(True)
	plt.xlabel(r'$' + xlabel + '$')
	plt.ylabel(r'$' + ylabel + '$')
	plt.title(run_name + '\n' + title)
	ax.set_xlim(float(min(times)), float(max(times)))

	i = 0
	while(i < len(timeseries_array)):
		ax.plot(times, timeseries_array[i], label=r'$' + label_array[i] + '$')
        	i+=1

	ax.legend(loc=2, prop=legendfont, borderaxespad=1.)

	full_fn = output_dir + '/' + filename
	plt.savefig(full_fn)

	print "Saved as " + full_fn

	plt.close()


# Compile all series to plot
i = 0
while i < ntrajcs:
	run = mylib.get_run_name()
	times = mylib.get_timeseries(i, 't')

	# All series to be plotted
	all_plot_series = {}

	# Get basic timeseries
	for srs in base_series:
		slabel = srs[0]
		sdata = mylib.get_timeseries(i, srs[1])

                all_plot_series[slabel] = sdata

	# Get summed timeseries
	for srs in summed_series:
		slabel = srs[0]

		slst = []
		for slbl in srs[1]:
			slst.append(all_plot_series[slbl])

		sdata = mylib.get_summed_timeseries(slst)

                all_plot_series[slabel] = sdata

	# Get derivative timeseries
        for srs in differentiated_series:
                slabel = srs[0]
                sdata = mylib.get_derivative_timeseries(all_plot_series[srs[1]])

                all_plot_series[slabel] = sdata

	# Get integrated timeseries
        for srs in integrated_series:
                slabel = srs[0]

		init_series = mylib.get_timeseries(i, srs[2])
		sdata = mylib.get_integrated_timeseries(all_plot_series[srs[1]], init_series[0])

		all_plot_series[slabel] = sdata
		
	# Get difference timeseries
        for srs in difference_series:
                slabel = srs[0]
                sdata = mylib.get_difference_timeseries(all_plot_series[srs[1][0]], all_plot_series[srs[1][1]])

                all_plot_series[slabel] = sdata

        # Get relative error timeseries
        for srs in relerror_series:
                slabel = srs[0]
                sdata = mylib.get_relerror_timeseries(all_plot_series[srs[1][0]], all_plot_series[srs[1][1]])

                all_plot_series[slabel] = sdata

	# Produce each requested plot
	for plot in plots:
                data = []
                labels = []

                # Get each quantity
                for srs in plot[0]:
                        qdata = all_plot_series[srs]
                        qlabel = srs

                        data.append(qdata)
                        labels.append(qlabel)

                units = plot[1]
                ylabel_base = plot[2]
                title_base = plot[3]
                fn_base = plot[4]

                title = 'trajc ' + str(i).zfill(5) + ' | ' + r'$' + title_base + '$'
                fn = fn_base + '_' + str(i).zfill(5) + '.png'
                xlabel = 't \ \ (s)'
                ylabel = ylabel_base + '\ \ ( ' + units + ')'

                calcplot( fn, times, data, labels, title, xlabel, ylabel )

	i+=1
