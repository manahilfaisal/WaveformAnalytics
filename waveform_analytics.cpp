#include <TGClient.h>
#include <TGFrame.h>
#include <TRootEmbeddedCanvas.h>
#include <TGSlider.h>
#include <TGLabel.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TRandom3.h>
#include <TMath.h>

#include <vector>
#include <algorithm>


void updateNoise(Int_t pos);


const int N = 1000;
double fs = 1000.0;   // Sampling frequency
double freq = 5.0;    // Sine frequency

TRandom3 rng(0);

std::vector<double> cleanSig(N), gaussSig(N), noisySig(N), corrSig(N);

TGraph *gClean = nullptr;
TGraph *gGauss = nullptr;
TGraph *gNoisy = nullptr;
TGraph *gCorr  = nullptr;

TGLabel *statsLabel  = nullptr;
TGLabel *sliderLabel = nullptr;

TRootEmbeddedCanvas *ec = nullptr;

double mean(const std::vector<double>& x) {
    double s = 0; for(auto v:x) s+=v; return s/x.size();
}

double median(std::vector<double> x) {
    std::sort(x.begin(),x.end());
    int n=x.size();
    return (n%2==0)?0.5*(x[n/2-1]+x[n/2]):x[n/2];
}

double mode(const std::vector<double>& x,int bins=50){
    double minV=*std::min_element(x.begin(),x.end());
    double maxV=*std::max_element(x.begin(),x.end());
    double bw=(maxV-minV)/bins;
    if(bw==0) return minV;

    std::vector<int> hist(bins,0);
    for(auto v:x){
        int idx=std::min(int((v-minV)/bw),bins-1);
        hist[idx]++;
    }
    int maxIdx=std::distance(hist.begin(), std::max_element(hist.begin(),hist.end()));
    return minV+bw*(maxIdx+0.5);
}

void crossCorrelation(){
    for(int lag=0;lag<N;lag++){
        double num=0, den1=0, den2=0;
        for(int i=0;i<N-lag;i++){
            num+=cleanSig[i]*noisySig[i+lag];
            den1+=cleanSig[i]*cleanSig[i];
            den2+=noisySig[i+lag]*noisySig[i+lag];
        }
        if(den1>0 && den2>0) corrSig[lag]=num/sqrt(den1*den2);
        else corrSig[lag]=0;
    }
}

void updateNoise(Int_t pos){
    double sigma = pos/100.0;
    if(sliderLabel)
        sliderLabel->SetText(TString::Format("Adjust Gaussian Noise σ = %.2f", sigma));

    // Update Gaussian & noisy signals
    for(int i=0;i<N;i++){
        gaussSig[i] = rng.Gaus(0,sigma);
        noisySig[i] = cleanSig[i] + gaussSig[i];

        gGauss->SetPoint(i,i/fs,gaussSig[i]);
        gNoisy->SetPoint(i,i/fs,noisySig[i]);
    }

    crossCorrelation();
    for(int i=0;i<N;i++)
        gCorr->SetPoint(i,i/fs,corrSig[i]);

    // Update statistics
    double mc=mean(cleanSig), mg=mean(gaussSig), mn=mean(noisySig);
    double medc=median(cleanSig), medg=median(gaussSig), medn=median(noisySig);
    double modec=mode(cleanSig), modeg=mode(gaussSig), moden=mode(noisySig);
    double corrPeak=*std::max_element(corrSig.begin(),corrSig.end());

    if(statsLabel)
        statsLabel->SetText(TString::Format(
            "Clean: μ=%.3f | median=%.3f | mode=%.3f\n"
            "Noise: μ=%.3f | median=%.3f | mode=%.3f\n"
            "Noisy: μ=%.3f | median=%.3f | mode=%.3f\n"
            "Correlation peak = %.3f",
            mc,medc,modec, mg,medg,modeg, mn,medn,moden, corrPeak
        ));

    // Redraw all updated graphs
    TCanvas *c = ec->GetCanvas();

    c->cd(2); gGauss->Draw("AL"); 
    c->cd(3); gNoisy->Draw("AL");
    c->cd(4); gCorr->Draw("AL");

    c->Modified();
    c->Update();
}

void noise_correlation(){

    // Generate clean sine
    for(int i=0;i<N;i++){
        double t=i/fs;
        cleanSig[i]=sin(2*TMath::Pi()*freq*t);
        gaussSig[i]=0;
        noisySig[i]=cleanSig[i];
        corrSig[i]=0;
    }
    crossCorrelation();

    TGMainFrame *f = new TGMainFrame(gClient->GetRoot(),1200,1000);

    ec = new TRootEmbeddedCanvas("ec",f,1200,800);
    f->AddFrame(ec,new TGLayoutHints(kLHintsExpandX|kLHintsExpandY));

    TCanvas *c = ec->GetCanvas();
    c->Divide(2,2); // 2x2 grid: Clean | Gaussian | Noisy | Correlation

    gClean = new TGraph(N);
    for(int i=0;i<N;i++) gClean->SetPoint(i,i/fs,cleanSig[i]);
    c->cd(1);
    gClean->SetLineColor(kBlue);
    gClean->SetTitle("Clean Signal;Time (s);Amplitude");
    gClean->Draw("AL");


    gGauss = new TGraph(N);
    for(int i=0;i<N;i++) gGauss->SetPoint(i,i/fs,gaussSig[i]);
    c->cd(2);
    gGauss->SetLineColor(kGreen+2);
    gGauss->SetTitle("Gaussian Noise;Time (s);Amplitude");
    gGauss->Draw("AL");


    gNoisy = new TGraph(N);
    for(int i=0;i<N;i++) gNoisy->SetPoint(i,i/fs,noisySig[i]);
    c->cd(3);
    gNoisy->SetLineColor(kRed);
    gNoisy->SetTitle("Noisy Signal;Time (s);Amplitude");
    gNoisy->Draw("AL");

    gCorr = new TGraph(N);
    for(int i=0;i<N;i++) gCorr->SetPoint(i,i/fs,corrSig[i]);
    c->cd(4);
    gCorr->SetLineColor(kBlue);
    gCorr->SetTitle("Normalized Cross-Correlation;Lag (s);Correlation");
    gCorr->Draw("AL");

    c->Update();


    TGHorizontalFrame *controls = new TGHorizontalFrame(f);

    sliderLabel = new TGLabel(controls,"Adjust Gaussian Noise σ = 0");
    controls->AddFrame(sliderLabel,new TGLayoutHints(kLHintsLeft,10,10,5,5));

    statsLabel = new TGLabel(controls,"Statistics will appear here");
    controls->AddFrame(statsLabel,new TGLayoutHints(kLHintsCenterX,10,10,5,5));

    TGHSlider *slider = new TGHSlider(controls,300,kSlider1|kScaleDownRight);
    slider->SetRange(0,300);
    slider->SetPosition(0);
    slider->Connect("PositionChanged(Int_t)",0,0,"updateNoise(Int_t)");
    controls->AddFrame(slider,new TGLayoutHints(kLHintsRight,10,10,5,5));

    f->AddFrame(controls,new TGLayoutHints(kLHintsExpandX));

    f->SetWindowName("Clean | Gaussian | Noisy | Correlation");
    f->MapSubwindows();
    f->Resize();
    f->MapWindow();
}
