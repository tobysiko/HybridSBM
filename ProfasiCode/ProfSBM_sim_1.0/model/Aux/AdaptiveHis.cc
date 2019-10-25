/*******************************************************************************
    PROFASI: Protein Folding and Aggregation Simulator, Version 1.5
    Copyright (C) (2012)  Anders Irback and Sandipan Mohanty
    Email: profasi@thep.lu.se
    Home Page: http://cbbp.thep.lu.se/activities/profasi/
    Version control (git) : https://trac.version.fz-juelich.de/PROFASI

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License
    (see PROFASI/gpl.txt).

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
********************************************************************************/

#include "AdaptiveHis.hh"
#include "profasi_io.hh"
#include <cstdio>
#include <math.h>
#include <unistd.h>

using prf::cerr;
using std::list;
using std::make_pair;
using std::pair;
using std::binary_function;
using std::vector;
using std::string;

using namespace prf;

namespace prf_utils
{
    AdaptiveHis::AdaptiveHis() : His1D(), adj_enabled(true) {}

    AdaptiveHis::AdaptiveHis(const AdaptiveHis &hs) : His1D(hs)
    {
        oor=hs.oor;
        adj_enabled=hs.adj_enabled;
    }

    AdaptiveHis & AdaptiveHis::operator=(const AdaptiveHis &hs)
    {
        if (this!=&hs) {
            His1D::operator=(hs);
            oor=hs.oor;
            adj_enabled=hs.adj_enabled;
        }

        return *this;
    }

    AdaptiveHis::~AdaptiveHis() {}

    void AdaptiveHis::init()
    {
        His1D::init();
        oor.clear();
    }

    void AdaptiveHis::reset()
    {
        His1D::reset();
        oor.clear();
    }

    int AdaptiveHis::put(double x, int i)
    {
        int ans=His1D::put(x,i);

        if (ans==0 && adj_enabled) oor.push_back(make_pair<double, int>(x,i));

        return ans;
    }

    int AdaptiveHis::nput(double howmanytimes, double x, int i)
    {
        int ans=His1D::nput(howmanytimes, x,i);

        if (ans==0 && adj_enabled) {
            for (int j=0;j<howmanytimes;++j)
                oor.push_back(make_pair<double, int>(x,i));
        }

        return ans;
    }

    int AdaptiveHis::adjust()
    {
        if (!adj_enabled) return 0;

        int i1=0,i2=nbin-1;

        find_new_bin_range(i1,i2);

        return adapt_bins(i1,i2);
    }

    int AdaptiveHis::adapt_bins(int minusedbin, int maxusedbin)
    {
        if (minusedbin!=0 or maxusedbin!=(nbin-1)) {
            xmin=xmin+minusedbin*xbin;
            xmax=xmax-(nbin-1-maxusedbin)*xbin;
            Matrix<double> oldxy=xy;
            int newnbin=maxusedbin-minusedbin+1;
            xy.allocate(nblk,newnbin);

            for (int i=0;i<nblk;++i) {
                for (int j=std::max(0,minusedbin);j<=std::min(maxusedbin,nbin-1);++j) {
                    xy[i][j-minusedbin]=oldxy[i][j];
                }
            }

            nbin=newnbin;

            list<pair<double, int> >::iterator it,jt;

            for (it=oor.begin(); it!=oor.end();) {
                int ibin=(int) floor((it->first - xmin)/xbin);
                jt=it; ++jt;

                if (ibin>=0 && ibin<newnbin) {
                    xy[it->second][ibin]++;
                    nentries_in[it->second]+=1;
                    oor.erase(it);
                }

                it=jt;
            }

            Logger(10)<<"AdaptiveHis["<<Name()<<"]> New range ("<<xmin<<", "<<xmax<<") with "
                <<nbin<<" bins. Num out of range = "<<oor.size()<<"\n";
            return 1;
        }

        return 0;
    }

    struct Comp :public binary_function<pair<double,int>,pair<double,int>,bool> {
        bool operator()(pair<double,int> &x, pair<double,int> &y) {
            return (x.first<y.first);
        }
    };

    int AdaptiveHis::find_new_bin_range(int &newmin, int &newmax)
    {
        Logger blog(20);
        int minusedbin=0, maxusedbin=(nbin-1);
        bool binunused=true;

        while (binunused && minusedbin<nbin) {
            for (int i=0;i<nblk;++i) {
                if (xy[i][minusedbin]!=0) {
                    binunused=false;
                    break;
                }
            }
            if (binunused) ++minusedbin; else break;
        }

        binunused=true;

        while (binunused && maxusedbin>=0) {
            for (int i=0;i<nblk;++i) {
                if (xy[i][maxusedbin]!=0) {
                    binunused=false;
                    break;
                }
            }
            if (binunused) --maxusedbin; else break;
        }

        int bmin(minusedbin),bmax(maxusedbin);

        if (!oor.empty()) sort_oorpoints();

        if (minusedbin==nbin || maxusedbin==-1) {
            prf::cerr<<"AdaptiveHis::find_new_bin_range()> None of the bins for "
                <<Name()<<" is populated. Current range is "<<xmin<<" to "<<xmax
                <<" with "<<nbin<<" bins\n";
            bmin=0; bmax=nbin-1;

            if (oor.empty()) {
                prf::cerr<<"The out of range histogram is also empty. "
                <<"Can not adjust range in this situation!\n";
                return 0;
            } else {
                list<pair<double,int> >::iterator it;
                double oormin,oormax;
                oormin=oormax=(oor.begin()->first);

                for (it=oor.begin(); it!=oor.end();++it) {
                    if (it->first < oormin) oormin=it->first;

                    if (it->first > oormax) oormax=it->first;
                }

                minusedbin=bmin=(int) floor((oormin-xmin)/xbin);

                maxusedbin=bmax=(int) floor((oormax-xmin)/xbin);
            }
        } else if (!oor.empty()) {
            blog<<"AdaptiveHis["<<Name()<<"]> Out of range list has "<<oor.size()<<" entries.\n";
            int nabove=0,nbelow=0;
            list<pair<double,int> >::iterator it;

            for (it=oor.begin(); it!=oor.end();++it) {
                if (it->first < xmin) ++nbelow;

                if (it->first >=xmax) ++nabove;
            }

            it=oor.begin();

            do {
                bmin=(int) floor((it->first-xmin)/xbin);

                if (nbelow<-bmin) {
                    ++it;--nbelow; bmin=minusedbin;
                } else break;
            } while (it!=oor.end());

            list<pair<double,int> >::reverse_iterator jt;

            jt=oor.rbegin();

            do {
                bmax=(int) floor((jt->first-xmin)/xbin);

                if (nabove<(bmax-(nbin-1))) {
                    ++jt;--nabove; bmax=maxusedbin;
                } else break;
            } while (jt!=oor.rend());
        }

        newmin=std::min(bmin,minusedbin);

        newmax=std::max(maxusedbin,bmax);

        if (newmin!=0 || newmax!=(nbin-1)) {
            // with a little bit of padding
            // reducing expensive calls to adjust
            newmin-=2;newmax+=2;
        }

        return 1;
    }

    int AdaptiveHis::sort_oorpoints()
    {
        if (!oor.empty()) {
            oor.sort(Comp());
        }

        /*    prf::cout<<"Sorted list of out of range points...\n";
            for (list<pair<double,int> >::iterator it=oor.begin();
            it!=oor.end(); ++it) {
                prf::cout<<it->first<<", "<<it->second<<"\n";
            }*/
        return 1;
    }

    AdaptiveHis add(vector<AdaptiveHis> &hs)
    {
        int nerr=0;

        if (hs.empty()) {
            prf::cerr<<"add: Empty list of AdaptiveHis\n";
            return AdaptiveHis();
        }

        double xmax=hs[0].Xmax(),xmin=hs[0].Xmin(),xbin=hs[0].Xbin();

        for (size_t i=1;i<hs.size();++i) nerr+=hs[0].check_compatibility(hs[i]);

        AdaptiveHis ans;

        if (nerr==0) {

            for (size_t i=0;i<hs.size();++i) {
                xmin=std::min(xmin,hs[i].Xmin())-xbin;
                xmax=std::max(xmax,hs[i].Xmax())+xbin;
            }

            int nbins=(int) (0.5+(xmax-xmin)/xbin);

            ans.Range(xmin,xmax);
            ans.Nbins(nbins);
            ans.NBlocks(hs[0].NBlocks());
            ans.init();

            for (size_t i=0;i<hs.size();++i) {
                for (int j=0;j<ans.NBlocks();++j) {
                    for (int k=0;k<hs[i].Nbins();++k) {
                        ans.nput(hs[i].data()[j][k],hs[i].xval(k),j);
                    }
                }
            }

            list<pair<double, int> >::iterator it;

            for (size_t i=0;i<hs.size();++i) {
                for (it=hs[i].out_of_range_list().begin();
                     it!=hs[i].out_of_range_list().end();++it) {
                    ans.put(it->first,it->second);
                }
            }

            ans.adjust();
        } else {
            prf::cerr<<"Addition of histograms failed because of "<<nerr
                    <<" errors\n";
        }

        return ans;
    }

    int AdaptiveHis::check_compatibility(const AdaptiveHis &hs)
    {
        int nerr=0;
        if (NBlocks()!=hs.NBlocks()) {
            prf::cerr<<"AdaptiveHis> "<<Name()<<": "<<NBlocks()<<" blocks and "
                    <<hs.Name()<<": "<<hs.NBlocks()<<" blocks.\n";
            ++nerr;
        }

        if (fabs((xbin-hs.Xbin())/(xbin))>0.0001) {
            prf::cerr<<"AdaptiveHis> Bin sizes differ. "<<Name()<<": "<<xbin
                    <<" and "<<hs.Name()<<": "<<hs.xbin<<"\n";
            ++nerr;
        }

        if (!is_good_bin_boundary(hs.Xmin())) {
            prf::cerr<<"AdaptiveHis> xmin value of "<<hs.Name()<<"("<<hs.Xmin()
                    <<") is not aligned at bin boundaries of "<<Name()
                    <<" with xmin = "<<xmin<<" and bin size = "<<xbin<<"\n";
            ++nerr;
        }
        if (!is_good_bin_boundary(hs.Xmax())) {
            prf::cerr<<"AdaptiveHis> xmax value of "<<hs.Name()<<"("<<hs.Xmax()
                    <<") is not aligned at bin boundaries of "<<Name()
                    <<" with xmin = "<<xmin<<" and bin size = "<<xbin<<"\n";
            ++nerr;
        }

        return nerr;
    }

    AdaptiveHis &AdaptiveHis::operator+=(AdaptiveHis &hs)
    {
        AdaptiveHis ans;
        ans.NBlocks(NBlocks());
        int nerr=check_compatibility(hs);
        if (nerr==0) {
            ans.Range(std::min(xmin,hs.Xmin())-xbin,std::max(xmax,hs.Xmax())+xbin);
            ans.Nbins((int)(((ans.Xmax()-ans.Xmin())/xbin)+0.5));
            ans.init();
            for (int j=0;j<ans.NBlocks();++j) {
                for (int k=0;k<Nbins();++k) {
                    ans.nput(data()[j][k],xval(k),j);
                }
                for (int k=0;k<hs.Nbins();++k) {
                    ans.nput(hs.data()[j][k],hs.xval(k),j);
                }
            }

            list<pair<double, int> >::iterator it;

            for (it=out_of_range_list().begin();
            it!=out_of_range_list().end();++it) {
                ans.put(it->first,it->second);
            }

            for (it=hs.out_of_range_list().begin();
            it!=hs.out_of_range_list().end();++it) {
                ans.put(it->first,it->second);
            }

            ans.adjust();
            (*this)=ans;
        }
        return *this;
    }

    AdaptiveHis add(vector<AdaptiveHis*> hsv)
    {
        int nerr=0;
        vector<AdaptiveHis*> hs;

        for (size_t i =0;i<hsv.size();++i) {
            if (hsv[i]!=NULL) hs.push_back(hsv[i]);
        }

        if (hs.empty()) {
            prf::cerr<<"add: Empty list of usable AdaptiveHis pointers\n";
            return AdaptiveHis();
        }

        double xmax=hs[0]->Xmax(),xmin=hs[0]->Xmin(),xbin=hs[0]->Xbin();

        for (size_t i=1;i<hs.size();++i) nerr+=hs[0]->check_compatibility(*hs[i]);

        AdaptiveHis ans;

        if (nerr==0) {

            for (size_t i=0;i<hs.size();++i) {
                xmin=std::min(xmin,hs[i]->Xmin())-xbin;
                xmax=std::max(xmax,hs[i]->Xmax())+xbin;
            }

            int nbins=(int) (0.5+(xmax-xmin)/xbin);

            ans.Range(xmin,xmax);
            ans.Nbins(nbins);
            ans.NBlocks(hs[0]->NBlocks());
            ans.init();

            for (size_t i=0;i<hs.size();++i) {
                for (int j=0;j<ans.NBlocks();++j) {
                    for (int k=0;k<hs[i]->Nbins();++k) {
                        ans.nput(hs[i]->data()[j][k],hs[i]->xval(k),j);;
                    }
                }
            }

            list<pair<double, int> >::iterator it;

            for (size_t i=0;i<hs.size();++i) {
                for (it=hs[i]->out_of_range_list().begin();
                     it!=hs[i]->out_of_range_list().end();++it) {
                    ans.put(it->first,it->second);
                }
            }

            ans.adjust();
        } else {
            prf::cerr<<"addition failed because of "<<nerr<<" errors\n";
        }

        return ans;
    }

    void AdaptiveHis::Export(const char *filename, int normmode, int lyout)
    {
        His1D::Export(filename,normmode,lyout);

        if (!oor.empty()) {
            Output fp((filename+string(".out_of_range")).c_str(),"w");

            for (list<pair<double, int> >::iterator it=oor.begin();
                 it!=oor.end();++it) {
                fp<<it->first<<"  "<<it->second<<"\n";
            }

            fp.close();
        } else unlink((filename+string(".out_of_range")).c_str());
    }

    bool AdaptiveHis::is_good_bin_boundary(double somept)
    {
        double eps=0.001;
        somept-=xmin;
        somept/=xbin;
        somept=fabs(somept);
        somept-=((int) somept);
        return (somept<eps or somept>(1-eps));
    }
}
