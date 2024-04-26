########################### PATH #############################

export PATH="$HOME/packages/postdoc/scripts:$PATH"
export PATH="$HOME/packages/rtRosetta/scripts:$PATH"
export PATH="$HOME/.gem/ruby/2.5.0/bin:$PATH"
export PATH="$HOME/.bin:$PATH" # .local/bin got polluted with random python crap
export PATH="$PATH:$HOME/packages/silent_tools"
export PATH="$PATH:$HOME/bii/swallow/scripts"

export PYTHONPATH="$HOME/packages:$PYTHONPATH"
export PYTHONPATH="$HOME/packages/NatureProtocols:$PYTHONPATH"
export PYTHONPATH="$HOME/packages/codon_harmony:$PYTHONPATH"
export PYTHONPATH="$HOME/packages/ppi_pipeline_tools:$PYTHONPATH"

######################### ALIASES #############################

alias l='ls -l --all --human-readable --no-group --color=auto --classify -v'
alias rsin='rsync -Rr --progress --update'
alias hist="sort | uniq -c | sort -r"
alias less='less -S'
alias tree='tree -C --filelimit=20'
alias rgf='rg --files | rg'
alias kol="column -s, -t"
alias watch='watch ' # triggers alias expansion
alias ptree='ps --user $(id -u) f'
alias ttime='/usr/bin/time -v'
alias countfiles='du -a | cut -d/ -f2 | sort | uniq -c | sort -nr'

######################### TERMINAL ############################

PS1='\[\e[38;5;178m\]\h\[\e[38;5;178m\]:\[\e[38;5;202m\]\w\[\e[38;5;99m\] ► \[\e[0m\]'

export TERM=xterm-color
set bell-style none
# can output a warning that messes with sftp
bind 'TAB: menu-complete' >/dev/null 2>/dev/null 

export EDITOR=vim
export VISUAL=$EDITOR

########################### DIGS #############################

alias sstatus='clusterstatus | (head -n 31; grep $USER)'
alias wstatus="watch 'clusterstatus | (head -n 31; grep $USER)'"
alias sq='squeue --user `whoami`'
alias sqh='sjob | less' 
alias fname='readlink -f'

if [ "$(hostname)" == "jojo" ]; then
    ulimit -u 10000
fi

source ~/bii/cloud/nextflow.sh

######################## COMPLETION ##########################

source ~/.fire_completion # python fire completion
# source /software/mmseqs2/util/bash-completion.sh

########################### EXTRA ############################

case "$-" in
*i*) # interactive shell
	bind 'set mark-symlinked-directories on'
	;;
*) # not interactive shell
	;;
esac

# load z for fuzzy cd from scripts/external/z.sh (this is scripts/.bashrc)

SCRIPTPATH="$( cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 ; pwd -P )"
Z=$SCRIPTPATH/external/z.sh
[ -f $Z ] && . $Z

weather(){ curl wttr.in/seattle_wa; }
cityonahill(){ curl wttr.in/san_francisco_ca; }
mouthofhell(){ curl wttr.in/cambridge_ma; }

nbless() {
    ipython nbconvert --to markdown --stdout $1 | mdless
}


# log bash input and search with hh
# nicked from wyang12's .bashrc
export PROMPT_COMMAND='echo "$(date "+%Y-%m-%d.%H:%M:%S") $(pwd) $(history 1)" >> ~/.logs_bash/$(date "+%Y-%m-%d").log;'
function hh() {
            ls -1 ~/.logs_bash/* | xargs grep -a $1 | tail -n ${2:-15}
}

export REMOTE=/home/wyang12/Documents/Binders/CTLA4/CTLA4_hits/L1_H1-3/2c_split_variants/final_split/2_split/

# csvkit uses tabulate, which is slow and lacks a streaming option
function csvless() {
    # csvlook=/home/dfeldman/.conda/envs/df-pyr-tf/bin/csvlook
    cat <(head -400 $1 | csvlook "${@:2}") <(csvlook $1 "${@:2}" | tail -n +402) | /usr/bin/less -S
}

# get typical digs paths
function real() {
    readlink -f $@ | sed 's/\/mnt//'
}

function wip() {
    cd `real ~/.wip`
}

function makewip() {
    if [ ! -z $1 ] 
    then 
        target=`real $1`
    else
        target=`real .`
    fi
    rm ~/.wip 2>/dev/null
    ln -s $target ~/.wip
}

# SO 11456403, allows unquoted wildcards like
# r path/to/stuff/X/
reset_expansion(){ CMD="$1";shift;$CMD "$@";set +f;}

function r {
    remote_path="$1"; shift
    rsync -rR --progress -z $@ ${DIGS_NODE:-jojo}:"${remote_path}" $HOME
}
alias r='set -f;reset_expansion r'

export BASH_SILENCE_DEPRECATION_WARNING=1
function rmnot {
    pat="$1"; shift
    ls -1 $@ | grep -v $pat | xargs rm
}


export MPLCONFIGDIR="$HOME/.config/matplotlib"


# sets title of terminal tab in jupyter lab
function title_terminal {
   PROMPT_COMMAND="echo -ne \"\033]0;$1\007\""
}

