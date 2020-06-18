export EDITOR=vim
export VISUAL=$EDITOR

SCRIPTPATH="$( cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 ; pwd -P )"
Z=$SCRIPTPATH/external/z.sh
[ -f $Z ] && . $Z

alias   l='ls -lg -h -a -G --color=auto'
alias   ls='ls -CF'
alias   rsync='rsync -Rr --progress --update'
alias   rename='rename -s'

PS1='\[\033[0;33m\]\u@\h \[\033[1;31m\]\w\[\033[${?/[^0]/39}m\]\$ \[\033[0;38m'

alias kol="column -s, -t"

set bell-style none
# can output a warning that messes with sftp
bind 'TAB: menu-complete' >/dev/null 2>/dev/null 

alias hist="sort | uniq -c | sort -r"
alias tar='tar -zxvf'
alias less='less -S'
alias rgf='rg --files | rg'

export TERM=xterm-color

export PYTHONPATH="$HOME/packages:$PYTHONPATH"
export PATH="$HOME/packages/postdoc/scripts:$PATH"
export PATH="$HOME/.gem/ruby/2.5.0/bin:$PATH"

weather(){ curl wttr.in/seattle_wa; }
cityonahill(){ curl wttr.in/san_francisco_ca; }
mouthofhell(){ curl wttr.in/cambridge_ma; }

alias sstatus='clusterstatus  | (head -n 33; grep $USER)'
alias wstatus="watch 'clusterstatus  | (head -n 33; grep $USER)'"
alias sq='squeue --user `whoami`'

nbless() {
    ipython nbconvert --to markdown --stdout $1 | mdless
}