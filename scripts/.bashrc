########################### PATH #############################

export PYTHONPATH="$HOME/packages:$PYTHONPATH"
export PATH="$HOME/packages/postdoc/scripts:$PATH"
export PATH="$HOME/.gem/ruby/2.5.0/bin:$PATH"

######################### ALIASES #############################

alias   l='ls -l --all --no-group --color=auto --classify -v'
alias   rsyncr='rsync -Rr --progress --update'
alias hist="sort | uniq -c | sort -r"
alias tar='tar -zxvf'
alias less='less -S'
alias rgf='rg --files | rg'
alias kol="column -s, -t"
# alias   rename='rename -s'

######################### TERMINAL ############################

PS1='\[\033[0;33m\]\u@\h \[\033[1;31m\]\w\[\033[${?/[^0]/39}m\]\$ \[\033[0;38m'
export TERM=xterm-color
set bell-style none
# can output a warning that messes with sftp
bind 'TAB: menu-complete' >/dev/null 2>/dev/null 

export EDITOR=vim
export VISUAL=$EDITOR

########################### DIGS #############################

alias sstatus='clusterstatus  | (head -n 31; grep $USER)'
alias wstatus="watch 'clusterstatus  | (head -n 31; grep $USER)'"
alias sq='squeue --user `whoami`'

########################### EXTRA ############################

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

