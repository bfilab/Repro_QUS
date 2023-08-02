function Uhk = UHK_loop(gama,alpha)

if round(alpha)==alpha
    Uhk1  = UHK( gama,alpha-1e-7 );
    Uhk2  = UHK( gama,alpha+1e-7 );
    Uhk = Uhk1 + (Uhk2-Uhk1)/2;
else
    Uhk = UHK(gama,alpha);
end

end