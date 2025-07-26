function S_2d=s4dto2d(S_4d)

for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
              if i==1 && j==1
               m=1;
                 elseif i==1 && j==2
                  m=6;
                 elseif i==1 && j==3
                  m=5;
                 elseif i==2 && j==1
                  m=6;
                 elseif i==2 && j==2
                  m=2;
                 elseif i==2 && j==3
                  m=4;
                 elseif i==3 && j==1
                  m=5;
                 elseif i==3 && j==2
                  m=4;
                 elseif i==3 && j==3
                  m=3;
              end
              if k==1 && l==1
               n=1;
                 elseif k==1 && l==2
                  n=6;
                 elseif k==1 && l==3
                  n=5;
                 elseif k==2 && l==1
                  n=6;
                 elseif k==2 && l==2
                  n=2;
                 elseif k==2 && l==3
                  n=4;
                 elseif k==3 && l==1
                  n=5;
                 elseif k==3 && l==2
                  n=4;
                 elseif k==3 && l==3
                  n=3;
              end
%% 1st equation
              if m==1 && n==1
                  S_2d(m,n)=S_4d(i,j,k,l);
              elseif m==1 && n==2
                  S_2d(m,n)=S_4d(i,j,k,l);
              elseif m==1 && n==3
                  S_2d(m,n)=S_4d(i,j,k,l);
              elseif m==2 && n==1
                  S_2d(m,n)=S_4d(i,j,k,l);
              elseif m==2 && n==2
                  S_2d(m,n)=S_4d(i,j,k,l);
              elseif m==2 && n==3
                  S_2d(m,n)=S_4d(i,j,k,l);
              elseif m==3 && n==1
                  S_2d(m,n)=S_4d(i,j,k,l);
              elseif m==3 && n==2
                  S_2d(m,n)=S_4d(i,j,k,l);
              elseif m==3 && n==3
                  S_2d(m,n)=S_4d(i,j,k,l);  
% 2nd cofficient 
              elseif m==1 && n==4
                  S_2d(m,n)=2*S_4d(i,j,k,l);
              elseif m==1 && n==5
                  S_2d(m,n)=2*S_4d(i,j,k,l);
              elseif m==1 && n==6
                  S_2d(m,n)=2*S_4d(i,j,k,l);
              elseif m==2 && n==4
                  S_2d(m,n)=2*S_4d(i,j,k,l);
              elseif m==2 && n==5
                  S_2d(m,n)=2*S_4d(i,j,k,l);
              elseif m==2 && n==6
                  S_2d(m,n)=2*S_4d(i,j,k,l);
              elseif m==3 && n==4
                  S_2d(m,n)=2*S_4d(i,j,k,l);
              elseif m==3 && n==5
                  S_2d(m,n)=2*S_4d(i,j,k,l);
              elseif m==3 && n==6
                  S_2d(m,n)=2*S_4d(i,j,k,l);
              elseif m==4 && n==1
                  S_2d(m,n)=2*S_4d(i,j,k,l);
              elseif m==4 && n==2
                  S_2d(m,n)=2*S_4d(i,j,k,l);
              elseif m==4 && n==3
                  S_2d(m,n)=2*S_4d(i,j,k,l);
              elseif m==5 && n==1
                  S_2d(m,n)=2*S_4d(i,j,k,l);
              elseif m==5 && n==2
                  S_2d(m,n)=2*S_4d(i,j,k,l);
              elseif m==5 && n==3
                  S_2d(m,n)=2*S_4d(i,j,k,l);
              elseif m==6 && n==1
                  S_2d(m,n)=2*S_4d(i,j,k,l);
              elseif m==6 && n==2
                  S_2d(m,n)=2*S_4d(i,j,k,l);
              elseif m==6 && n==3
                  S_2d(m,n)=2*S_4d(i,j,k,l);
% 4nd cofficient 
              elseif m==4 && n==4
                  S_2d(m,n)=4*S_4d(i,j,k,l);
              elseif m==4 && n==5
                  S_2d(m,n)=4*S_4d(i,j,k,l);
              elseif m==4 && n==6
                  S_2d(m,n)=4*S_4d(i,j,k,l);
              elseif m==5 && n==4
                  S_2d(m,n)=4*S_4d(i,j,k,l);
              elseif m==5 && n==5
                  S_2d(m,n)=4*S_4d(i,j,k,l);
              elseif m==5 && n==6
                  S_2d(m,n)=4*S_4d(i,j,k,l);
              elseif m==6 && n==4
                  S_2d(m,n)=4*S_4d(i,j,k,l);
              elseif m==6 && n==5
                  S_2d(m,n)=4*S_4d(i,j,k,l);
              elseif m==6 && n==6
                  S_2d(m,n)=4*S_4d(i,j,k,l);
              end
           
            end
        end
     end 
 end
end