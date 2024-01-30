%================================Compute and save the binomial PDF table
function prep_binomial_table(r_max,N_end,min_variant_read,min_total_read,SFS_totalsteps_base,option_ploidy)
%-----------------------------------------Prepare the binomial PDF table
matrix_binomial_PDF = zeros(r_max,N_end,SFS_totalsteps_base);
for r=1:r_max
    fprintf('%d/%d\n',r,r_max);
    for m=1:N_end
        for i=1:SFS_totalsteps_base
            if (r<min_variant_read)||(r<min_total_read)
                matrix_binomial_PDF(r,m,i)  = 0;
                continue;
            end
            %               Find boundaries for s: s is in ( (i-1)*r/SFS_totalsteps_base,i*r/SFS_totalsteps_base ]
            r1  = r*(i-1)/SFS_totalsteps_base;
            if r1<ceil(r1)
                r1  = ceil(r1);
            else
                r1  = r1+1;
            end
            r2  = floor(r*i/SFS_totalsteps_base);
            %               Find Prob(s/r is in ((i-1)/SFS_totalsteps_base,i/SFS_totalsteps_base] given m and r)
            Prob= 0;
            for s=max(min_variant_read,r1):r2
                Prob    = Prob+binopdf(s,r,m/(N_end*option_ploidy));
            end
            matrix_binomial_PDF(r,m,i)  = Prob;
        end
    end
end
%--------------------------------------------Save the binomial PDF table
filename=['Binomial_PDF_' num2str(N_end) ...
    '_' num2str(r_max) ...
    '_' num2str(min_variant_read) ...
    '_' num2str(min_total_read) ...
    '_' num2str(SFS_totalsteps_base) ...
    '_' num2str(option_ploidy) ...
    '.mat'];
save(filename,'matrix_binomial_PDF');
end