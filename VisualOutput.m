function VisualOutput(model,Diagnostics)

    tW= model.t2/model.dt/model.do;
    sC= get(groot,'Screensize');
    
    tX= ceil( [ 0,0.25,0.5,0.75,1 ]*tW );
    lX= ceil( [ 0,0.25,0.5,0.75,1 ]*model.t2/model.dt );
    
    tE= [ 1,ceil( [ 0.25,0.5,0.75 ]*model.Ne ),model.Ne+1 ];
    lE= [ 1,ceil( [ 0.25,0.5,0.75 ]*model.Ne ),model.Ne+1 ];

    bL= [ 51,153,255 ]/255;
    rD= [ 255,51,51 ]/255;
    gR= [ 0,204,0 ]/255;
    
    set( 0,'DefaulttextFontName','Arial' )
    
    
    %%% STATE %%%
    %%% ===== %%%
    figure('uni','pi','pos',[ sC(3)*0.05,sC(4),sC(3)*0.55,sC(4)*0.40 ])
    subplot(231)
    plot( Diagnostics(1).mseXf,'-b' ); hold on 
    plot( Diagnostics(1).espXf,'-r' ); grid on
    plot( Diagnostics(1).eskXf,'-g' )
    
    bMi= max( 0,Diagnostics(1).espXf(end)-mean(Diagnostics(1).espXf) ); 
    bMa= round( Diagnostics(1).espXf(1)+mean(Diagnostics(1).espXf),1 );
    
    set(gca,'FontSize',15,'XLim',[0,tW],'XTick',tX,'XTickLabel',lX,'YLim', ...
        [ bMi,bMa ],'YTick',linspace( bMi,bMa,4 ),'YTickLabelRotation',90 )
    
    title('State variables: Experiment #1','FontSize',18)
    xlabel('Modeling steps','FontSize',16); ylabel('Forcast ensemble statistics','FontSize',16);
    legend('Mean-Squared-Errors','Ensemble-Spread','Prediction-Skill','Location','Best')
    
    subplot(232)
    x= plot( Diagnostics(1).dXf,'-b' ); hold on 
    y= plot( mean( Diagnostics(1).dXf,2 ),'-.r','LineWidth',1.2 ); grid on
    z= plot( Diagnostics(1).refX,'--k','LineWidth',1.2 );
    
    set(gca,'FontSize',15,'XLim',[tW*0.25,tW*0.5],'XTick',tX,'XTickLabel',lX )
    
    title( [ 'State variable #', num2str(model.Ix) ],'FontSize',18)
    xlabel('Modeling steps','FontSize',16); ylabel('Forcast ensemble','FontSize',16);
    legend([ x(1),y,z ],'Forecast-Ensemble-Members','Forecast-Ensemble-Mean','Reference-Trajectory','Location','Best')
    
    subplot(233)
    plot(Diagnostics(1).skwXf,'Color',bL); hold on 
    plot(Diagnostics(1).kurXf,'Color',rD); grid on
    
    set(gca,'FontSize',15,'XLim',[0,tW],'XTick',tX,'XTickLabel',lX )
    
    title( 'High Order Moments','FontSize',18)
    xlabel('Modeling steps','FontSize',16); ylabel('Forcast ensemble','FontSize',16);
    legend('Skewness','Kurtosis','Location','Best')
    
    subplot(234)
    [ ax,s1h1,s1h2 ] = plotyy( 1:tW,Diagnostics(1).inbXf,1:tW,Diagnostics(1).indXf );
    axes(ax(1)); hold on; inb= plot( 1:tW,ones( tW,1 ),':b','LineWidth',1 ) ;
    axes(ax(2)); hold on; ind= plot( 1:tW,sqrt(model.Ne+1)/sqrt(2*model.Ne)*ones( tW,1 ),':r','LineWidth',1 ) ;
    set( s1h1,'LineStyle','-','Color',bL,'LineWidth',1.2 ); 
    set( s1h2,'LineStyle','-','Color',rD,'LineWidth',1.2 );
    set( ax(1),'YColor','k','FontSize',15,'XGrid','on','XLim',[0,tW],'XTick',tX,'XTickLabel',lX,'YTickLabelRotation',90 )
    set( ax(2),'YColor','k','FontSize',15,'XLim',[0,tW],'XTick',[],'YLim',[ 0.5,ceil(Diagnostics(1).indXf(1)) ], ...
               'YTick',linspace( 0.5,ceil(Diagnostics(1).indXf(1)),length(get(ax(1),'YTick')) ),'YTickLabelRotation',90 )
    set( get(ax(1),'Ylabel'),'String','Filter inbreeding criterion','FontSize',16 )
    set( get(ax(1),'Xlabel'),'String','Modeling Steps','FontSize',16 )
    set( get(ax(2),'Ylabel'),'String','Consistency ratio','FontSize',16 )
    title( 'Ensemble Reliability Measures','FontSize',18 )
    legend( [ s1h1,s1h2,inb,ind ],'Inbreeding','Consistency','Ideal Performance','Ideal Performance','Location','Best' )
    
    subplot(235)
    H= rankHist( Diagnostics(1).dXf,Diagnostics(1).refX );
    bar( (1:model.Ne+1),H*100 ); grid on
    
    set(gca,'FontSize',15,'XLim',[1,model.Ne+1],'XTick',tE,'XTickLabel',lE,'YLim',[0,10] );
    title( [ 'Rank-Histogram for state variable #', num2str(model.Ix) ],'FontSize',18)
    xlabel('bins','FontSize',16); ylabel('Frequency (%)','FontSize',16);
    
    subplot(236)
    semilogx( 1:tW,Diagnostics(1).dfsXf,'Color',bL,'LineWidth',1.2 ); hold on
    semilogx( 1:tW,1-Diagnostics(1).dfsXf,'Color',rD,'LineWidth',1.2 ); grid on
    
    set(gca,'FontSize',15,'XLim',[0,tW],'XTick',[ 1,10,100,tW(end) ],'XTickLabel',[ 1,10,100,tW(end) ]*model.do )
    title( 'Analysis Sensitivity','FontSize',18)
    xlabel('Modeling steps','FontSize',16); ylabel('Forcast ensemble','FontSize',16);
    legend( 'Observation Influence','Background Influence','Location','Best' )
    
    
    
    %%% PARAMETERS %%%
    %%% ========== %%%
    figure('uni','pi','pos',[ sC(3)*0.05,sC(4)*0.05,sC(3)*0.55,sC(4)*0.40 ])
    subplot(231)
    plot( Diagnostics(1).mseFf,'-b' ); hold on 
    plot( Diagnostics(1).espFf,'-r' ); grid on
    plot( Diagnostics(1).eskFf,'-g' )
    
    bMi= max( 0,Diagnostics(1).espFf(end)-mean(Diagnostics(1).espFf) ); 
    bMa= round( Diagnostics(1).espFf(1)+mean(Diagnostics(1).espFf),1 );
    
    set(gca,'FontSize',15,'XLim',[0,tW],'XTick',tX,'XTickLabel',lX,'YLim', ...
        [ bMi,bMa ],'YTick',linspace( bMi,bMa,4 ),'YTickLabelRotation',90 )
    
    title('Parameters: Experiment #1','FontSize',18)
    xlabel('Modeling steps','FontSize',16); ylabel('Forcast ensemble statistics','FontSize',16);
    legend('Mean-Squared-Errors','Ensemble-Spread','Prediction-Skill','Location','Best')
    
    subplot(232)
    x= plot( Diagnostics(1).dFf,'-b' ); hold on 
    y= plot( mean( Diagnostics(1).dFf,2 ),'-.r','LineWidth',1.2 ); grid on
    z= plot( Diagnostics(1).refF,'--k','LineWidth',1.2 );
    
    set(gca,'FontSize',15,'XLim',[0,tW],'XTick',tX,'XTickLabel',lX )
    
    title( [ 'Parameter #', num2str(model.Ip) ],'FontSize',18)
    xlabel('Modeling steps','FontSize',16); ylabel('Forcast ensemble','FontSize',16);
    legend([ x(1),y,z ],'Forecast-Ensemble-Members','Forecast-Ensemble-Mean','Reference-Parameter-Value','Location','Best')
    
    subplot(233)
    plot(Diagnostics(1).skwFf,'Color',bL); hold on 
    plot(Diagnostics(1).kurFf,'Color',rD); grid on
    
    set(gca,'FontSize',15,'XLim',[0,tW],'XTick',tX,'XTickLabel',lX )
    
    title( 'High Order Moments','FontSize',18)
    xlabel('Modeling steps','FontSize',16); ylabel('Forcast ensemble','FontSize',16);
    legend('Skewness','Kurtosis','Location','Best')

    subplot(234)
    kobs= ceil( [ 1/tW,0.01,0.05,0.1,1 ]*tW );
    prob= zeros( length(kobs),100 );
    xbin= zeros( length(kobs),100 );
    for i= 1:length(kobs)
        [ prob(i,:),xbin(i,:) ]= ksdensity( Diagnostics(1).dFf(kobs(i),:) );
    end
    plot( xbin(1,:),prob(1,:),'-k' ); hold on
    plot( xbin(2,:),prob(2,:),'-b' ); grid on
    plot( xbin(3,:),prob(3,:),'-r' ); 
    plot( xbin(4,:),prob(4,:),'-g' ); 
    plot( xbin(5,:),prob(5,:),'-m' ); aa= axis; 
    plot( Diagnostics(1).refF(1).*ones(1,100),linspace(aa(3),aa(4),100),'--k' );
    
    set(gca,'FontSize',15 );
    title( [ 'Pobability Density Function, parameter #', num2str(model.Ip) ],'FontSize',18)
    xlabel('Parameter Value','FontSize',16); ylabel('PDF','FontSize',16);
    legend( [ 'Step #' num2str( model.do*kobs(1) ) ],[ 'Step #' num2str( model.do*kobs(2) ) ], ...
            [ 'Step #' num2str( model.do*kobs(3) ) ],[ 'Step #' num2str( model.do*kobs(4) ) ], ...
            [ 'Step #' num2str( model.do*kobs(5) ) ],'Truth','Location','Best' )
        
    subplot(235)    
    cdf= ksdensity( Diagnostics(1).dFf(kobs(i),:),'function','cdf' ); hold on
    z1=plot( xbin(i,:),cdf,'-b' ); grid on
    z2=plot( mean(xbin(i,:))*ones(1,100),cdf,'-.r' ); aa= axis; box on
    z3=plot( xbin(i,80)*ones(1,100),linspace(aa(3),aa(4),100),'-.','Color',gR);
    plot( xbin(i,20)*ones(1,100),linspace(aa(3),aa(4),100),'-.','Color',gR);
    z5=plot( Diagnostics(1).refF(1).*ones(1,100),linspace(aa(3),aa(4),100),'--k' );
    
    set(gca,'FontSize',15 );
    title( [ 'Cumulative Distribution, parameter #', num2str(model.Ip) ],'FontSize',18)
    xlabel('Parameter Value','FontSize',16); ylabel('PDF','FontSize',16);
    legend( [ z1,z2,z3,z5 ],[ 'Step #' num2str( model.do*kobs(5) ) ],'Sample-Mean','80% Bound','Truth','Location','Best' )
    
    subplot(236)
    stairs( 1:tW,Diagnostics(1).corFf,'Color',gR,'LineWidth',1.2 ); grid on
    
    set(gca,'FontSize',15,'XLim',[0,tW],'XScale','log','XTick',[ 1,10,100,tW(end) ],'XTickLabel',[ 1,10,100,tW(end) ]*model.do )
    
    title( 'Cross-correlations with Observed Variables','FontSize',18)
    xlabel('Modeling steps','FontSize',16); ylabel('Normalized forecast statistics','FontSize',16);
    
    
    %%% TEXT OUTPUT %%%
    %%% =========== %%%
    disp( [ '> STATE: Average of "Mean-Squared-Errors" = ' num2str( round( mean(Diagnostics(1).mseXf),3 ) ) ] )
    disp( [ '> STATE: Average of "Ensemble-Spread"     = ' num2str( round( mean(Diagnostics(1).espXf),3 ) ) ] )
    disp( [ '> STATE: Average of "Prediction-Skill"    = ' num2str( round( mean(Diagnostics(1).eskXf),3 ) ) ] )
    
    disp( ' ' )
    
    disp( [ '> PARAMETERS: Average of "Mean-Squared-Errors" = ' num2str( round( mean(Diagnostics(1).mseFf),3 ) ) ] )
    disp( [ '> PARAMETERS: Average of "Ensemble-Spread"     = ' num2str( round( mean(Diagnostics(1).espFf),3 ) ) ] )
    disp( [ '> PARAMETERS: Average of "Prediction-Skill"    = ' num2str( round( mean(Diagnostics(1).eskFf),3 ) ) ] )
    
    
    %%% SAVE OUTPUT %%%
    %%% =========== %%%
    save( [ model.na '_N' num2str( model.Ne ) '_p' num2str( model.do ) ...
            '_E' num2str( model.Ns )  '_' model.aa '_' model.us '.mat' ],'model','Diagnostics')
    
end