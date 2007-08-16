      subroutine equaly(n1,n2,obs1,year1,obs2,year2,obsday1,obsday2,
     +count)
      
      implicit none
      
      integer n1,n2,obs1(*),obs2(*),year1(*),year2(*),count,x,y
      integer obsday1(*),obsday2(*)
      
      count = 0
      
      do x=1,n1
      	 do y=1,n2
            if(year1(x).eq.year2(y)) then
	       count = count + 1
	       obsday1(count) = obs1(x)
	       obsday2(count) = obs2(y)
	    endif
         enddo
      enddo
      
      return
      end
