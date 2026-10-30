for i=4:length(Fobj)
    A1(i)=( abs(Fobj(i))-abs(Fobj(i-1)) )/(abs(Fobj(i-1)));
    A2(i)=( abs(Fobj(i))-abs(Fref) )/(abs(Fref));
    A3(i)=abs(min(A2(i),0));
    fi(i)=A1(i)+5*A3(i);
end