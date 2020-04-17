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
bind 'TAB: menu-complete'

alias hist="sort | uniq -c | sort -r"
alias tar='tar -zxvf'
alias less='less -S'
alias rgf='rg --files | rg'

export TERM=xterm-color

export PYTHONPATH="$HOME/packages:$PYTHONPATH"

weather(){ curl wttr.in/seattle_wa; }
cityonahill(){ curl wttr.in/san_francisco_ca; }
mouthofhell(){ curl wttr.in/cambridge_ma; }
