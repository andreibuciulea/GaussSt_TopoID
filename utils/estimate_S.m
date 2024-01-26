function [S_hat, out] = estimate_S(C,model,reg)
%reg.S_true= mbinarize(C-diag(diag(C)),2);
    switch model
        case 'gti-GST'
            out = gti_GST(C,reg);
            S_hat = out.S;
        case 'GSR'
            [S_hat, out]= GSR(C,reg);
        case 'GGSR'
            [S_hat, out]= GGSR(C,reg);
        case 'GGSR-fast'
            [S_hat, out]= GGSR_fast(C,reg);
        case 'GGSR-bi'
            out = gti_bi_GST(C,reg);
            S_hat = out.S;
        case 'GGSR-j'
            out = GGSR_j(C,reg);
            S_hat = out.S_hat;
        case 'GGSR-j-v1'
            out = GGSR_j_v1(C,reg);
            S_hat = out.S_hat;
        case 'GSR-noise'
            [S_hat, out]= GSR_noise(C,reg);
        case 'rec-adj'
            [S_hat, out]= rec_adj(C,reg);
        case 'est-S-G'
            [S_hat, out]= est_S_G(C,reg);

        case 'est-S-GST-sparse'
            [S_hat, out]= est_S_GST_sparse(C,reg);
        case 'est-S-GST-block'
            [S_hat, out]= est_S_GST_block(C,reg);
        case 'est-S-GST-3B'
            [S_hat, out]= est_S_GST_3blck(C,reg);
        case 'est-S-GST-3BI'
            [S_hat, out]= est_S_GST_3BI(C,reg);
        case 'est-S-GST-3BRw'
            [S_hat, out]= est_S_GST_3BRw(C,reg);
        case 'est-S-GST-fast'
            [S_hat, out]= est_S_GST_fast(C,reg);

        case 'est-S-GST-3BFast'
            [S_hat, out]= est_S_GST_3BFast(C,reg);
        case 'est-S-GST-3B-noise'
            [S_hat, out]= est_S_GST_3B_noise(C,reg);
        case 'est-S-GST-noise'
            [S_hat, out]= est_S_GST_noise(C,reg);
        case 'GL'
            [S_hat, out] = graphicalLasso(C,reg);
            S_hat = S_hat - diag(diag(S_hat)); 
        otherwise
            error(['Error: unknown estimation method' model])
    end
end
