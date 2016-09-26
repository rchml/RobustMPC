
    spaceDim = 4;
    N = 10;
    Lnum = 1;
    Lmin = .2;
    Lmax = 3.4;
    Cnum = 1;
    Cmin = .8;
    Cmax = 6.2;
    betaPol = 1;
    gamPol = 1;

    s1 = sprintf('DATA_SIM/exp_xR%d_N%d_L%d[%0.1f,%0.1f]_C%d[%0.1f,%0.1f]_BetaPol%d_GamPol%d', ...
            spaceDim, N, Lnum, Lmin, Lmax, Cnum, Cmin, Cmax, betaPol, gamPol);

    if size(dir(s1),1) == 0
        mkdir(s1);
    else
        fprintf('Warning: writing experiments in an old directory\n');
        pause
    end
    
    %%
    
    Lx = Lmin + .3;
    Cy = Cmin + .8;
    
    s2 = strcat(s1, sprintf('/L%0.1f_C%0.1f', Lx, Cy));
    
    if size(dir(s2),1) == 0
        mkdir(s2);
    else
        fprintf('Warning: writing experiments in an old directory\n');
        pause
    end
    
    %%
    
    saveas(gcf,strcat(s1, sprintf('/xdata.png')));
    savefig(gcf, strcat(s1, sprintf('/xdata')));
    
    saveas(gcf,strcat(s1, sprintf('/Betadata.png')));
    
    %%
    
    a.N = N;
    a.s1 = s1;
    save(strcat(s1, sprintf('/simdata.mat')), '-struct', 'a');
    %%
    
    string = s1;
    fid = fopen(strcat(s1, sprintf('/output.txt')), 'wt');
    fprintf(fid, string);
    fclose(fid);
    
    