function x = DhDNtot(h,lambda,n)

    a = chiExact(h,lambda,n);
    b = dmudlam(h,lambda,n);
    c = dchidh(h,lambda,n);
    d = dchidlam(h,lambda,n);
    e = dmudN(h, lambda, n);
    f = dchidN(h, lambda, n);
    x = -((-d *e + b *f)/(b *c - a* d));

end