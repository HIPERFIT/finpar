#!/usr/bin/env python
import numpy as N
import pylab as P
import sys

_CONTRACT_LIST=[1,3,2]
axes_ylabels=["GPU SpeedUp: gaming vs. mobile solution"] # ["SpeedUp with nVidia GTX 680 GPU","SpeedUp with nVidia Quadro 2000M GPU"]
_DATA_SCALINGS=[1.0,1.0]

_BARS_TEXT_FONTSIZES=[14,12] # 10
_BARS_TEXT_OFFSETS=[20,12]
_HATCHES=[' ',"x"]
_AXES_YLIM_INCREASE_RATIO = 1.115
_AXES_XTICKS_FONTSIZE=16
_AXES_YLABEL_FONTSIZE=16
_TITLE_FONTSIZE=18
_TEXT_EXTRA_FONTSIZE=15
_TEXT_LEGEND_FONTSIZE=_TEXT_EXTRA_FONTSIZE
_FIGURE_SIZE=(8.85,6.6) # inches
_FIGURE_DPI=96

_FIGURE_SAVE=True
_GRAYSCALE=False
_FIGURE_EXT=".png" # ".eps"/".pdf"/".png"

# .eps? no transparency. Use gray scale
if _GRAYSCALE==True or _FIGURE_EXT=='.eps':
  _TRANSPARENCY_ALPHAS=[1.0]*4
  _BARS_COLORS=['0.6','0.9','0.75','0.4']
  _TEXT_COLORS=["k"]*2
else:
  _TEXT_COLORS=["b","r"]
  _TRANSPARENCY_ALPHAS=[0.9,0.9,0.6,0.6]
  _BARS_COLORS=['b','y','r','c']
#------------------------------------------------
_DEBUG=5
#=========================================================================
def timing_file_read(filename):
  # dims:
  # opt_set, contract, measures (2d)
  # dict, dict, list of tuples of 2 elemens
  with open(filename,"r") as fH:
    timings_read=eval(fH.read())
  return timings_read
#---------------------------------------------------------------
def timings_parse(timings):
  db={}
  for opt_set in sorted(timings.keys()):
    # read timings per contract
    db[opt_set]={}
    for contract in sorted(timings[opt_set].keys()):
      measure=timings[opt_set][contract]
      if not len(measure): measure=[(1.0,0.0)]
      measure=N.array(measure)
      # TODO: point estimate is fine, std() needs better estimation
      measure=measure[:,1]/measure[:,0]
      measure_avg=measure.mean()
      measure_std=measure.std()
      db[opt_set][contract]=[measure_avg,measure_std]
      if _DEBUG>=10:
        print "# DEBUG: opt_set: %d, contract: %d -> %8.3f +/- %7.3f" % (
          opt_set,contract,measure_avg,measure_std
          )
  return db
#-------------------------------------------------------------------------
def flatten(l):
  "flattens a list l"
  import operator
  return reduce(operator.add, l)
#-------------------------------------------------------------------------
def image_sub_ax(
    db,
    ax,
    opt_set_list,  # benchmarking configurations
    contract_list, # list of contracts to show
    width_total,   # width of a group of bars
    text_color='b',
    tick_color='b',
    error_color='k',
    axis_label_color='k',
    bar_posLeft_offset=0,
    bars_text_offset=20,
    bar_text_fontsize=12,
    bars_colors=_BARS_COLORS,
    transparency_alphas=_TRANSPARENCY_ALPHAS,
    scaling=1.0,     # Note: this scaling should be used in development only
    hatch='',
    autolabel_flag=False,
    ylabel="ylabel",
    text_rotation=90,
    yticks_visible=False,
    ax_ylim=True,
    ):
  shapes_opt_set=[]
  max_val=0
  for opt_set in opt_set_list:
    bar_width=width_total/len(opt_set)
    shapes_opt_set=[None]*len(opt_set)
    for (contract_j,contract) in enumerate(contract_list):
      for opt_j,opt in enumerate(opt_set):
        avg=db[opt][contract][0]*scaling
        std=db[opt][contract][1]*scaling
        if _DEBUG>=50: print "opt: %d, contract: %d -> %.1f +/- %.1f" % (opt,contract,avg,std)  
        # update max value
        if max_val<avg: max_val=avg
        # bars positioning
        bar_left=bar_posLeft_offset+(contract_j+1)+opt_j*bar_width
        # return a bit left if any of the first two bars were zero
        if opt_j>=1:
          if db[opt_set[0]][contract][0]==0: bar_left-=bar_width/2
          if db[opt_set[1]][contract][0]==0: bar_left-=bar_width/2
        # bars
        bar_height=avg
        shape=ax.bar(
          left=bar_left,
          height=bar_height,
          width=bar_width,
          color=bars_colors[opt_j],
          yerr=std,
          ecolor=error_color,
          alpha=transparency_alphas[opt_j],
          hatch=hatch,
        )
        # legend
        if shapes_opt_set[opt_j]==None:
          shapes_opt_set[opt_j]=shape
        # text with values over bars
        if autolabel_flag and bar_height>0: # and (opt_j%2==0):
          ax.text(bar_left+bar_width/2., bar_height+bars_text_offset, '%d'%int(float(bar_height)/scaling), ha='center', va='bottom', color=text_color, fontsize=bar_text_fontsize,
            rotation=text_rotation, backgroundcolor="w")
    bar_posLeft_offset+=len(contract_list)
  # ytick text color
  for t in ax.get_yticklabels():
    if yticks_visible==True: t.set_color(tick_color)
    else: t.set_visible(False)
  # ax: ylabel
  ax.set_ylabel(ylabel,color=axis_label_color,fontsize=_AXES_YLABEL_FONTSIZE)
  # ax: ylim
  # ax.set_ylim(0,ax_ylim)
  if ax_ylim==True: ax.set_ylim(0,max_val*_AXES_YLIM_INCREASE_RATIO)
  #
  return shapes_opt_set
#-------------------------------------------------------------------------
def image_draw(
    db,
    contract_list=_CONTRACT_LIST,
    opt_set_list=[ [0,1,2,3],[4,5,6,7] ],
    width_total=0.9,
    xtick_labels=[['SimpleF', 'MediumF', 'ComplexF'],['SimpleD','MediumD','ComplexD']],
    legend_labels=['a','b','c','d','A','B','C','D'],
    axes_ylabels=["ax1","ax2"],
    title="title",
    figure_extra_text_flags=['gaming','mobile'],
    scalings=_DATA_SCALINGS,
    grid=False,
    bars_colors=_BARS_COLORS,
    transparency_alphas=_TRANSPARENCY_ALPHAS,
    bars_text_offsets=_BARS_TEXT_OFFSETS,
    text_colors=_TEXT_COLORS,
    hatches=_HATCHES,
    bar_posLeft_offsets=[0.0,0.0],
    text_rotations=[90,0],
    filename_pre=None,
    extra_text_pos=(0.7,0.7),
    bar_text_fontsizes=_BARS_TEXT_FONTSIZES,
    axis_label_color='k',
  ):
  fig=P.figure(figsize=_FIGURE_SIZE,dpi=_FIGURE_DPI)
  # multiple axes?
  ax1=fig.add_subplot(111)
  ax_list=[ax1]
  if len(axes_ylabels)==2:
    ax2=ax1.twinx()
    ax_list.append(ax2)
  # bars for each ax, or collapsed in the same ax
  ax_j=0
  shapes=[]
  ylim_set=[True,False]
  for db_j,db_content in enumerate(db):
    # skip if this database is empty
    if db_content==None: continue
    #
    if len(axes_ylabels)>1:
      ax_j=db_j
    # bars
    shapes.append(
      image_sub_ax(
        db=db_content,
        ax=ax_list[ax_j],
        contract_list=contract_list,
        opt_set_list=opt_set_list,
        autolabel_flag=True,
        bars_colors=bars_colors,
        transparency_alphas=transparency_alphas,
        axis_label_color=axis_label_color,
        #
        scaling=scalings[db_j],
        bars_text_offset=bars_text_offsets[db_j],
        text_color=text_colors[db_j],
        hatch=hatches[db_j],
        text_rotation=text_rotations[db_j],
        ax_ylim=ylim_set[db_j],
        #
        width_total=width_total,
        ylabel=axes_ylabels[ax_j],
        tick_color=text_colors[ax_j],
        bar_text_fontsize=bar_text_fontsizes[ax_j],
        #
        bar_posLeft_offset=bar_posLeft_offsets[ax_j],
        )
    )
    # yticks: hide
    ytl = ax1.get_yticklabels()
    ytl[0].set_visible(False)
  # legend
  P.legend( shapes[0], legend_labels, prop={'size':_TEXT_LEGEND_FONTSIZE} )
  # xtick labels
  xticks_offset=width_total/2
  if len(axes_ylabels)==2: xticks_offset=width_total
  P.xticks(
    (N.arange(len(flatten(xtick_labels)))+xticks_offset+1),
    (flatten(xtick_labels)),
    fontsize=_AXES_XTICKS_FONTSIZE
    )
  # grid?
  if grid==True: ax1.grid(axis="y")
  # title
  P.title(title,fontsize=_TITLE_FONTSIZE)
  # figure extra text
  figure_extra_text(flags=figure_extra_text_flags,ax=ax1,extra_text_pos=extra_text_pos,axis_label_color=axis_label_color)
  # save?
  if _FIGURE_SAVE==True and filename_pre:
    filename=filename_pre+_FIGURE_EXT
    fig.savefig(filename,dpi=_FIGURE_DPI,bbox_inches='tight',pad_inches=0.1)
    if _DEBUG>=5: print "# DEBUG: saving '%s'" % filename

def figure_extra_text(flags=['gaming','mobile'],
    ax=None,
    extra_text_pos=(0.7,0.7),axis_label_color='k',
    ):
  if 'mobile' in flags:
    P.text(extra_text_pos[0]+0.005,extra_text_pos[1],       "XX",color=axis_label_color,transform = ax.transAxes,fontsize=_TEXT_EXTRA_FONTSIZE)
    P.text(extra_text_pos[0]+0.05, extra_text_pos[1]+0.003, ": nVidia Quadro 2000M",color=axis_label_color,transform = ax.transAxes,fontsize=_TEXT_EXTRA_FONTSIZE)
  if 'gaming' in flags:
    P.text(extra_text_pos[0]+0.015,extra_text_pos[1]+0.090, "XXX",color=axis_label_color,rotation='vertical',transform = ax.transAxes,fontsize=_TEXT_EXTRA_FONTSIZE)
    P.text(extra_text_pos[0]+0.05, extra_text_pos[1]+0.068, ": nVidia GeForce GTX 680",color=axis_label_color,transform = ax.transAxes,fontsize=_TEXT_EXTRA_FONTSIZE)
#=========================================================================
# MAIN
#=========================================================================
if __name__=="__main__":
  if len(sys.argv)<2:
    print "# USAGE: %s log_file1 [log_file2: optional]" % sys.argv[0]
    sys.exit(1)
  # first file log
  db=[]
  filename_list=sys.argv[1:3]
  for filename_j,filename in enumerate(filename_list):
    if _DEBUG>=5: print "# DEBUG: reading '%s'" % filename
    timings=timing_file_read(filename)
    # append timings, in command line order
    db.append( timings_parse(timings) )
  if _DEBUG>=100: print "# DEBUG: db\n%s" % (db,)
  # images: CG-Vect.
  image_draw(
    [db[0],None],
    contract_list=_CONTRACT_LIST,
    opt_set_list=[ [0,2],[4,6] ],
    legend_labels=['Coarse Grained All Opt ON','Vectorized All Opt ON']*2,
    axes_ylabels=["GPU SpeedUp: gaming solution"],
    figure_extra_text_flags=['gaming'],
    title="Switching Coarse-Grained/Vectorized Optimization\nSpeedup w.r.t. Sequential CPU Runtime (-O3)",
    scalings=_DATA_SCALINGS,
    bars_text_offsets=[8,2],
    filename_pre="optimizations-GPU-CG-0",
    extra_text_pos=(0.53,0.595),
    bars_colors=[_BARS_COLORS[0],_BARS_COLORS[2]],
    transparency_alphas=[_TRANSPARENCY_ALPHAS[0],_TRANSPARENCY_ALPHAS[2]],
    hatches=[' ']*2,
    )
  image_draw(
    [None,db[1]],
    contract_list=_CONTRACT_LIST,
    opt_set_list=[ [0,2],[4,6] ],
    legend_labels=['Coarse Grained All Opt ON','Vectorized All Opt ON']*2,
    axes_ylabels=["GPU SpeedUp: mobile solution"],
    figure_extra_text_flags=['mobile'],
    title="Switching Coarse-Grained/Vectorized Optimization\nSpeedup w.r.t. Sequential CPU Runtime (-O3)",
    scalings=_DATA_SCALINGS,
    bars_text_offsets=[8,2],
    filename_pre="optimizations-GPU-CG-1",
    extra_text_pos=(0.53,0.595),
    bars_colors=[_BARS_COLORS[0],_BARS_COLORS[2]],
    transparency_alphas=[_TRANSPARENCY_ALPHAS[0],_TRANSPARENCY_ALPHAS[2]],
    hatches=[' ']*2,
    )
  # images: SR
  image_draw(
    [db[0],None],
    contract_list=_CONTRACT_LIST,
    opt_set_list=[ [0,1,2,3],[4,5,6,7] ],
    legend_labels=['Coarse Grained All Opt ON','Coarse Grained SR OFF','Vectorized All Opt ON','Vectorized SR OFF']*2,
    axes_ylabels=["GPU SpeedUp: gaming solution"],
    figure_extra_text_flags=['gaming'],
    title="Switching ON/OFF Strength-Reduction (SR) Optimization\nSpeedup w.r.t. Sequential CPU Runtime (-O3)",
    scalings=_DATA_SCALINGS,
    bars_text_offsets=[8,2],
    filename_pre="optimizations-GPU-SR-0",
    extra_text_pos=(0.53,0.595),
    hatches=[' ']*2,
    )
  image_draw(
    [None,db[1]],
    contract_list=_CONTRACT_LIST,
    opt_set_list=[ [0,1,2,3],[4,5,6,7] ],
    legend_labels=['Coarse Grained All Opt ON','Coarse Grained SR OFF','Vectorized All Opt ON','Vectorized SR OFF']*2,
    axes_ylabels=["GPU SpeedUp: mobile solution"],
    figure_extra_text_flags=['mobile'],
    title="Switching ON/OFF Strength-Reduction (SR) Optimization\nSpeedup w.r.t. Sequential CPU Runtime (-O3)",
    scalings=_DATA_SCALINGS,
    bars_text_offsets=[8,2],
    filename_pre="optimizations-GPU-SR-1",
    extra_text_pos=(0.53,0.595),
    hatches=[' ']*2,
    )
  """
  # images: BD
  image_draw(
    db,
    contract_list=_CONTRACT_LIST,
    opt_set_list=[ [0,8,2,9],[4,10,6,11] ],
    legend_labels=['Coarse Grained All Opt ON','Coarse Grained BD OFF','Vectorized All Opt ON','Vectorized BD OFF']*2,
    axes_ylabels=axes_ylabels,
    #title="Switching ON/OFF Branch-Divergence (BD) Optimization\nThe Other Two Optimizations Are ON\nSpeedup w.r.t. Sequential CPU Runtime (-O3)",
    title="Switching ON/OFF Branch-Divergence (BD) Optimization\nSpeedup w.r.t. Sequential CPU Runtime (-O3)",
    scalings=_DATA_SCALINGS,
    filename_pre="optimizations-GPU-BD",
    extra_text_pos=(0.53,0.595),
    )
  """
  # images: MC
  image_draw(
    [db[0],None],
    contract_list=_CONTRACT_LIST,
    opt_set_list=[ [2,13],[6,15] ],
    legend_labels=['Vectorized All Opt ON','Vectorized MC OFF']*2,
    axes_ylabels=["GPU SpeedUp: gaming solution"],
    title="Switching ON/OFF Memory-Coalescence (MC) Optimization\nSpeedup w.r.t. Sequential CPU Runtime (-O3)",
    figure_extra_text_flags=['gaming'],
    scalings=_DATA_SCALINGS,
    bars_text_offsets=[8,1],
    filename_pre="optimizations-GPU-MC-0",
    extra_text_pos=(0.55,0.70),
    bars_colors=_BARS_COLORS[2:],
    transparency_alphas=_TRANSPARENCY_ALPHAS[2:],
    hatches=[' ']*2,
    )
  image_draw(
    [None,db[1]],
    contract_list=_CONTRACT_LIST,
    opt_set_list=[ [2,13],[6,15] ],
    legend_labels=['Vectorized All Opt ON','Vectorized MC OFF']*2,
    axes_ylabels=["GPU SpeedUp: mobile solution"],
    figure_extra_text_flags=['mobile'],
    title="Switching ON/OFF Memory-Coalescence (MC) Optimization\nSpeedup w.r.t. Sequential CPU Runtime (-O3)",
    scalings=_DATA_SCALINGS,
    bars_text_offsets=[8,1],
    filename_pre="optimizations-GPU-MC-1",
    extra_text_pos=(0.55,0.70),
    bars_colors=_BARS_COLORS[2:],
    transparency_alphas=_TRANSPARENCY_ALPHAS[2:],
    hatches=[' ']*2,
    )
  #
  P.show()
