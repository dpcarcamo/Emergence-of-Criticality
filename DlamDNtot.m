function x = DlamDNtot(h,lambda,n)

    %y -> -((-c e + a f)/(b c - a d))

    a = chiExact(h,lambda,n);
    b = dmudlam(h,lambda,n);
    c = dchidh(h,lambda,n);
    d = dchidlam(h,lambda,n);
    e = dmudN(h, lambda, n);
    f = dchidN(h, lambda, n);
    x = -((c* e - a* f)/(b* c - a *d));
end