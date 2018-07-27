unsigned int ViewToDAQChan(unsigned int ViewChan)
{
  size_t crate = ViewChan / 320;
  size_t Chan311;

  ViewChan = 8*(ViewChan/8+1)-ViewChan%8 -1;

  if(crate == 0)
  {
    ViewChan = 32*(ViewChan/32+1)-ViewChan%32 -1;
    size_t card = 4 - ((ViewChan / 32) % 5);
    if(ViewChan > 159)
    {
        size_t shift = 31 - (ViewChan % 32);
        Chan311 = (2*card)*32 + shift;
    }
    else
    {
       size_t shift = 31 - (ViewChan % 32);
       Chan311 = (2*card + 1)*32 + shift;
    }
  }
  else
  {
     size_t new_ViewChan = ViewChan - crate*320;
     size_t card = ((new_ViewChan / 32) % 5);
     if(new_ViewChan > 159)
     {
        size_t shift = new_ViewChan % 32;
        Chan311 = (2*card)*32 + shift;
     }
     else
     {
       size_t shift = new_ViewChan % 32;
       Chan311 = (2*card + 1)*32 + shift;
     }
     Chan311 = Chan311 + crate*320;
  } // end of if/else statementi

  return Chan311;
}
