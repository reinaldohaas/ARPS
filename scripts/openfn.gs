* Open multiple files
*
* @Author Yunheng Wang  11/07/2003
*

function main(args)

  say '------------- Opening file(s) -------------------'
  say ' *** ' args ' ***'

  count=1
  fn = subwrd(args, count)
  while(fn != "")
    'open ' fn
    count = count+1
    fn = subwrd(args, count)
  endwhile
  if (rc != 0)
    say '------------ Error while opening ----------------'
  else
    say '------------ ' count-1 ' files opened ---------------------'
  endif
return
