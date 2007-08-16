      subroutine combinations(id,n,ncomb,stat1,stat2)
      
      implicit none
      
      integer id(*),n,ncomb,x,y,count
      integer stat1(*),stat2(*)
      
      count = 0
      
      do x=1,n-1
      	 do y=x+1,n
            count = count + 1
            stat1(count) = id(x)
            stat2(count) = id(y)
         enddo
      enddo
      
      return
      end
