      subroutine altdep(x,y,alt,obsday,n,delBB)
    
      implicit none

      real pigrad, tot_weight, tot_dist_weight
      real argument, q_altobs, distance, dif_alt, dif_obs
      real dist_weight, local_weight
      integer istart, alt(*), obsday(*), n, i, j
      double precision x(*), y(*), delBB

      pigrad = 3.14159 / 180
      tot_weight = 0.
      tot_dist_weight = 0.
      
            
      do istart=1,n-1
        if(alt(istart).lt.100) then
	  goto 99
	endif
          do j=istart+1,n
             dif_alt = alt(istart) - alt(j)
              dif_obs = obsday(istart) - obsday(j)
                q_altobs = dif_obs / dif_alt
         argument = sin(x(istart) * pigrad) * (sin(x(j) * pigrad))+ 
     +		 cos(x(istart) * pigrad) * cos(x(j) * pigrad) * 
     +           cos(abs(y(istart) * pigrad) - (y(j) * pigrad))
		distance = 6371.11 * acos(argument)
               if(distance.lt.5) then
                 distance = 5
               endif
              dist_weight = 1.0 / distance
              if(alt(j).lt.100.or.abs(dif_alt).lt.50) then
                 dist_weight = 0.
              endif
             tot_dist_weight = tot_dist_weight + dist_weight
            local_weight = q_altobs * dist_weight
            if(alt(j).lt.100.or.abs(dif_alt).lt.50) then
               local_weight = 0.
            endif
           tot_weight = tot_weight + local_weight
          enddo
 99       continue
      enddo
      delBB = tot_weight / tot_dist_weight

      return
      end 
