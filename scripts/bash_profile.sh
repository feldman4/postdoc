export EDITOR=vim
export VISUAL=$EDITOR

alias   l='ls -lg -h -a -G'
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
alias rgf='rg --file | rg'

export TERM=xterm-color

export PYTHONPATH="$HOME/packages:$PYTHONPATH"

weather(){ curl wttr.in/seattle_wa; }
cityonahill(){ curl wttr.in/san_francisco_ca; }
mouthofhell(){ curl wttr.in/cambridge_ma; }
