* Open multiple HDF control files
*
* @Author Yunheng Wang  11/07/2003
*
function main(args)

  say '--------- Opening HDF control file(s) -------------'
  say ' *** ' args ' ***'

  count=1
  fn = subwrd(args, count)
  while(fn != "")
    'xdfopen ' fn
    count = count+1
    fn = subwrd(args, count)
  endwhile
  if (rc != 0)
    say '------------ Error while opening ----------------'
  else
    say '------------ ' count-1 ' files opened ---------------------'
  endif
  'set mproj off'
return
