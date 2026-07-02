function [u,x,dudl] = muchi(h,lambda,N)

    % u = real(muExact(h,lambda,N));
    % x = real(chiExact(h,lambda,N));
    % 
    % dudl = real(M12(h,lambda,N)/N);

    [u,x,dudl] = muChiExact2Spin(h, lambda, N);
end