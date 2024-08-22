function  Init_Env(fontsize, img_max_size)
    set(0,'DefaultFigureWindowStyle','docked')
    set(0,'defaultAxesFontSize',fontsize)
    set(groot, 'defaultLegendInterpreter','latex');
    set(groot, 'defaultAxesTickLabelInterpreter','latex');
    set(groot, 'defaulttextinterpreter','latex');  
    if img_max_size == 1
        set(groot, 'defaultFigureUnits','normalized')
        set(groot, 'defaultFigurePosition',[0 0 1 1])
    end
        
end

