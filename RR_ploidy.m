function rr=RR_ploidy(H_real,H_rec)
[R,L]=size(H_rec);
error=PloidMEC(H_rec,H_real);
rr=1-error/(L*R);
