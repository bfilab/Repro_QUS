function Xhk = XHK_loop(gama,alpha)

if round(alpha)==alpha
    Xhk1  = XHK( gama,alpha-1e-7 );
    Xhk2  = XHK( gama,alpha+1e-7 );
    Xhk = Xhk1 + (Xhk2-Xhk1)/2;
else
    Xhk = XHK(gama,alpha);
end

end