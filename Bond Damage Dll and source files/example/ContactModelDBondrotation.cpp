// contactmodeldbondrotation.cpp
#include "contactmodeldbondrotation.h"

#include "module/interface/icontactmechanical.h"
#include "module/interface/icontact.h"
#include "fish/src/parameter.h"
#include "module/interface/ipiecemechanical.h"
#include "module/interface/ipiece.h"
#include "C:\Program Files\Itasca\PFC700\PluginFiles\contactmodel\example/version1.txt"
#include "base/src/basetoqt.h"

#include "module/interface/ifishcalllist.h"
#include "utility/src/tptr.h"
#include "shared/src/mathutil.h"
#include "math.h"
#include "algorithm"
#include "stdlib.h"
#include "kernel/interface/iprogram.h"
#include "module/interface/icontactthermal.h"
#include "contactmodel/src/contactmodelthermal.h"




#ifdef DBondrotation_LIB
#ifdef _WIN32
  int __stdcall DllMain(void *,unsigned, void *)
  {
    return 1;
  }
#endif
  extern "C" EXPORT_TAG const char *getName()
  {
#if DIM==3
    return "contactmodelmechanical3ddbondrotation";
#else
    return "contactmodelmechanical2ddbondrotation";
#endif
  }

  extern "C" EXPORT_TAG unsigned getMajorVersion()
  {
    return MAJOR_VERSION;
  }

  extern "C" EXPORT_TAG unsigned getMinorVersion()
  {
    return MINOR_VERSION;
  }

  extern "C" EXPORT_TAG void *createInstance()
  {
    cmodelsxd::ContactModelDBondrotation *m = NEWC(cmodelsxd::ContactModelDBondrotation());
    return (void *)m;
  }
#endif // DBONDROTATION_EXPORTS

namespace cmodelsxd {
    static const quint32 linKnMask      = 0x00002; // Base 1!
    static const quint32 linKsMask      = 0x00004;
    static const quint32 linFricMask    = 0x00008;

    using namespace itasca;

    int ContactModelDBondrotation::index_ = -1;
    UInt ContactModelDBondrotation::getMinorVersion() const { return MINOR_VERSION;}

    ContactModelDBondrotation::ContactModelDBondrotation() : inheritanceField_(linKnMask|linKsMask|linFricMask)
                                                       , kn_(0.0)
                                                       , ks_(0.0)
                                                       , fric_(0.0)
		                                               , delta_lamda_(0.0)
	                                                   , prev_D_(0.0)
		                                               , prev_dellam_(0.0)
		                                               , normal_total_disp_(0.0) 
		                                               , shear_total_dispy_(0.0)
		                                               , shear_total_dispz_(0.0)
	                                                   , shear_total_disp_(0.0)
		                                               , bending_total_roy_(0.0)
		                                               , bending_total_roz_(0.0)
	                                                   , bending_total_ro_(0.0)                  
	                                                   , normal_plastic_disp_(0.0) 
	                                                   , shear_plastic_disp_(0.0)
	                                                   , bending_plastic_ro_(0.0)
		                                               , BendingMoy_(0.0)
		                                               , BendingMoz_(0.0)
		                                               , BendingMo_(0.0)
		                                               , Vy_(0.0)
		                                               , Vz_(0.0)
	                                                   , prev_npd_(0.0)                          
	                                                   , prev_spd_(0.0)
	                                                   , prev_bpr_(0.0) 
		                                               , npd1_(0.0)
		                                               , spd1_(0.0)
		                                               , bpr1_(0.0)
	                                                   , force_return_norm_(0.0)                 
	                                                   , force_return_tany_(0.0)                  
	                                                   , force_return_tanz_(0.0)                  
	                                                   , moment_return_by_(0.0)                   
	                                                   , moment_return_bz_(0.0)                    
                                                       , lin_F_(DVect(0.0))
                                                       , lin_S_(false)
                                                       , lin_mode_(0)
                                                       , rgap_(0.0)
                                                       , dpProps_(0)
                                                       , pbProps_(0)
                                                       , userArea_(0)
                                                       , energies_(0)
                                                       , effectiveTranslationalStiffness_(DVect2(0.0)) {
//    setFromParent(ContactModelMechanicalList::instance()->find(getName()));
    }

    ContactModelDBondrotation::~ContactModelDBondrotation() {
        if (dpProps_)
            delete dpProps_;
        if (pbProps_)
            delete pbProps_;
        if (energies_)
            delete energies_;
    }

    void ContactModelDBondrotation::archive(ArchiveStream &stream) {
        stream & kn_;
        stream & ks_;
        stream & fric_;
        stream & lin_F_;
        stream & lin_S_;
        stream & lin_mode_;
		stream & delta_lamda_;
		stream & prev_D_;
		stream & prev_dellam_;
		stream & normal_total_disp_;
		stream & shear_total_dispy_;
		stream & shear_total_dispz_;
		stream & shear_total_disp_;
		stream & bending_total_roy_;
		stream & bending_total_roz_;
		stream & bending_total_ro_;
		stream & normal_plastic_disp_;
		stream & shear_plastic_disp_;
		stream & bending_plastic_ro_;
		stream & BendingMoy_;
		stream & BendingMoz_;
		stream & BendingMo_;
		stream & Vy_;
		stream & Vz_;
		stream & prev_npd_;
		stream & prev_spd_;
		stream & prev_bpr_;
		stream & npd1_;
		stream & spd1_;
		stream & bpr1_;
		stream & force_return_norm_;
		stream & force_return_tany_;
		stream & force_return_tanz_;
		stream & moment_return_by_;
		stream & moment_return_bz_;

        if (stream.getArchiveState()==ArchiveStream::Save) {
            bool b = false;
            if (dpProps_) {
                b = true;
                stream & b;
                stream & dpProps_->dp_nratio_;
                stream & dpProps_->dp_sratio_;
                stream & dpProps_->dp_mode_;
                stream & dpProps_->dp_F_;
            }
            else
                stream & b;

            b = false;
            if (energies_) {
                b = true;
                stream & b;
                stream & energies_->estrain_;
                stream & energies_->eslip_;
                stream & energies_->edashpot_;
                stream & energies_->epbstrain_;
            }
            else
                stream & b;

            b = false;
            if (pbProps_) {
                b = true;
                stream & b;
                stream & pbProps_->pb_state_;
                stream & pbProps_->pb_rmul_;
                stream & pbProps_->pb_kn_;
                stream & pbProps_->pb_ks_;
				stream & pbProps_->pb_kn1_;
				stream & pbProps_->pb_ks1_;
                stream & pbProps_->pb_mcf_;
                stream & pbProps_->pb_ten_;
				stream & pbProps_->pb_com_;
				//stream & pbProps_->pb_mb_;
				stream & pbProps_->pb_ten1_;
				stream & pbProps_->pb_com1_;
				//stream & pbProps_->pb_mb1_;
                stream & pbProps_->pb_coh_;
				stream & pbProps_->pb_coh1_;
                //stream & pbProps_->pb_fa_;
                stream & pbProps_->pb_F_;
                stream & pbProps_->pb_M_;
				stream & pbProps_->mn_;
				stream & pbProps_->ms_;
				stream & pbProps_->mm_;
				stream & pbProps_->Fyield_;
				stream & pbProps_->ucn_;
				stream & pbProps_->ucn1_;
				stream & pbProps_->ucs_;
				stream & pbProps_->ucb_;
				//stream & pbProps_->alpha_;
				//stream& pbProps_->timing_;
				stream & pbProps_->D_;
				stream & pbProps_->Fy_;
				
            }
            else
                stream & b;

        } else {
            bool b(false);
            stream & b;
            if (b) {
                if (!dpProps_)
                    dpProps_ = NEWC(dpProps());
                stream & dpProps_->dp_nratio_;
                stream & dpProps_->dp_sratio_;
                stream & dpProps_->dp_mode_;
                stream & dpProps_->dp_F_;
            }
            stream & b;
            if (b) {
                if (!energies_)
                    energies_ = NEWC(Energies());
                stream & energies_->estrain_;
                stream & energies_->eslip_;
                stream & energies_->edashpot_;
                stream & energies_->epbstrain_;
            }
            stream & b;
            if (b) {
                if (!pbProps_)
                    pbProps_ = NEWC(pbProps());
				stream & pbProps_->pb_state_;
				stream & pbProps_->pb_rmul_;
				stream & pbProps_->pb_kn_;
				stream & pbProps_->pb_ks_;
				stream & pbProps_->pb_kn1_;
				stream & pbProps_->pb_ks1_;
				stream & pbProps_->pb_mcf_;
				stream & pbProps_->pb_ten_;
				stream & pbProps_->pb_com_;
				//stream & pbProps_->pb_mb_;
				stream & pbProps_->pb_ten1_;
				stream & pbProps_->pb_com1_;
				//stream & pbProps_->pb_mb1_;
				stream & pbProps_->pb_coh_;
				stream & pbProps_->pb_coh1_;
				//stream & pbProps_->pb_fa_;
				stream & pbProps_->pb_F_;
				stream & pbProps_->pb_M_;
				stream & pbProps_->mn_;
				stream & pbProps_->ms_;
				stream & pbProps_->mm_;
				stream & pbProps_->Fyield_;
				stream & pbProps_->ucn_;
				stream & pbProps_->ucn1_;
				stream & pbProps_->ucs_;
				stream & pbProps_->ucb_;
				//stream & pbProps_->alpha_;
				//stream & pbProps_->timing_;
				stream & pbProps_->D_;
				stream & pbProps_->Fy_;
				
            }
        }

        stream & inheritanceField_;
        stream & effectiveTranslationalStiffness_;

        if (stream.getArchiveState()==ArchiveStream::Save || stream.getRestoreVersion() == getMinorVersion())
            stream & rgap_;

        if (stream.getArchiveState() == ArchiveStream::Save || stream.getRestoreVersion() > 1)
            stream & userArea_;
    }

    void ContactModelDBondrotation::copy(const ContactModel *cm) {
        ContactModelMechanical::copy(cm);
        const ContactModelDBondrotation *in = dynamic_cast<const ContactModelDBondrotation*>(cm);
        if (!in) throw std::runtime_error("Internal error: contact model dynamic cast failed.");
        kn(in->kn());
        ks(in->ks());
        fric(in->fric());

		Dellam(in->Dellam());
		Prev_dvar(in->Prev_dvar());
		Prev_dellam(in->Prev_dellam());
		Normal_total_disp(in->Normal_total_disp());
		Shear_total_dispy(in->Shear_total_dispy());
		Shear_total_dispz(in->Shear_total_dispz());
		Shear_total_disp(in->Shear_total_disp());
		Bending_total_roy(in->Bending_total_roy());
		Bending_total_roz(in->Bending_total_roz());
		Bending_total_ro(in->Bending_total_ro());
		Normal_plastic_disp(in->Normal_plastic_disp());
		Shear_plastic_disp(in->Shear_plastic_disp());
		Bending_plastic_ro(in->Bending_plastic_ro());
		BendingMoy(in->BendingMoy());
		BendingMoz(in->BendingMoz());
		BendingMo(in->BendingMo());
		Vy(in->Vy());
		Vz(in->Vz());
		Prev_npd(in->Prev_npd());
		Prev_spd(in->Prev_spd());
		Prev_bpr(in->Prev_bpr());
		Npd1(in->Npd1());
		Spd1(in->Spd1());
		Bpr1(in->Bpr1());
		ReFnorm(in->ReFnorm());
		ReFtany(in->ReFtany());
		ReFtanz(in->ReFtanz());
		ReMby(in->ReMby());
		ReMbz(in->ReMbz());

        lin_F(in->lin_F());
        lin_S(in->lin_S());
        lin_mode(in->lin_mode());
        rgap(in->rgap());

        if (in->hasDamping()) {
            if (!dpProps_)
                dpProps_ = NEWC(dpProps());
            dp_nratio(in->dp_nratio());
            dp_sratio(in->dp_sratio());
            dp_mode(in->dp_mode());
            dp_F(in->dp_F());
        }
        if (in->hasEnergies()) {
            if (!energies_)
                energies_ = NEWC(Energies());
            estrain(in->estrain());
            eslip(in->eslip());
            edashpot(in->edashpot());
            epbstrain(in->epbstrain());
        }
        if (in->hasPBond()) {
            if (!pbProps_)
                pbProps_ = NEWC(pbProps());
            pbState(in->pbState());
            pbRmul(in->pbRmul());
            pbKn(in->pbKn());
            pbKs(in->pbKs());
			pbKn1(in->pbKn1());
			pbKs1(in->pbKs1());
            pbMCF(in->pbMCF());
            pbTen(in->pbTen());
			pbCom(in->pbCom());
			//pbMb(in->pbMb());
			pbTen1(in->pbTen1());
			pbCom1(in->pbCom1());
			//pbMb1(in->pbMb1());
            pbCoh(in->pbCoh());
			pbCoh1(in->pbCoh1());
            //pbFA(in->pbFA());
            pbF(in->pbF());
            pbM(in->pbM());
			Mn(in->Mn());
			Ms(in->Ms());
			Mm(in->Mm());
			Yield(in->Yield());
			Ucn(in->Ucn());
			Ucn1(in->Ucn1());
			Ucs(in->Ucs());
			Ucb(in->Ucb());
			//Alpha(in->Alpha());
			//Timing(in->Timing());
			Dvar(in->Dvar());
			Fyi(in->Fyi());
            pbTransStiff(in->pbTransStiff());
            pbAngStiff(in->pbAngStiff());
        }
        userArea_ = in->userArea_;
        inheritanceField(in->inheritanceField());
        effectiveTranslationalStiffness(in->effectiveTranslationalStiffness());
    }

    QVariant ContactModelDBondrotation::getProperty(uint i,const IContact *con) const {
        QVariant var;
        switch (i) {
        case kwLinKn:        return kn_;
        case kwLinKs:        return ks_;
        case kwLinFric:      return fric_;
        case kwLinMode:      return lin_mode_;
        case kwLinF:         var.setValue(lin_F_); return var;
        case kwLinS:         return lin_S_;
        case kwRGap:         return rgap_;
		case kwNtd:          return normal_total_disp_;
		case kwStdy:         return shear_total_dispy_;
		case kwStdz:         return shear_total_dispz_;
		case kwStd:          return shear_total_disp_;
		case kwBtry:         return bending_total_roy_;
		case kwBtrz:         return bending_total_roz_;
	    case kwBtr:          return bending_total_ro_;
		case kwNpd:          return normal_plastic_disp_;
		case kwSpd:          return shear_plastic_disp_;
		case kwBpr:          return bending_plastic_ro_;
		case kwBendingMy:    return BendingMoy_;
		case kwBendingMz:    return BendingMoz_;
		case kwBendingM:     return BendingMo_;
		case kwSheary:       return Vy_;
		case kwShearz:       return Vz_;
        case kwEmod: {
                        const IContactMechanical *c(convert_getcast<IContactMechanical>(con));
                        if (c ==nullptr) return 0.0;
                        double rsq(std::max(c->getEnd1Curvature().y(),c->getEnd2Curvature().y()));
                        double rsum(0.0);
                        if (c->getEnd1Curvature().y())
                            rsum += 1.0/c->getEnd1Curvature().y();
                        if (c->getEnd2Curvature().y())
                            rsum += 1.0/c->getEnd2Curvature().y();
                        if (userArea_) {
#ifdef THREED
                            rsq = std::sqrt(userArea_ / dPi);
#else
                            rsq = userArea_ / 2.0;
#endif
                            rsum = rsq + rsq;
                            rsq = 1. / rsq;
                        }
#ifdef TWOD
                        return (kn_ * rsum * rsq / 2.0);
#else
                        return (kn_ * rsum * rsq * rsq) / dPi;
#endif
                    }
        case kwKRatio:      return (ks_ == 0.0) ? 0.0 : (kn_/ks_);
        case kwDpNRatio:    return dpProps_ ? dpProps_->dp_nratio_ : 0;
        case kwDpSRatio:    return dpProps_ ? dpProps_->dp_sratio_ : 0;
        case kwDpMode:      return dpProps_ ? dpProps_->dp_mode_ : 0;
        case kwUserArea:    return userArea_;
        case kwDpF: {
                dpProps_ ? var.setValue(dpProps_->dp_F_) : var.setValue(DVect(0.0));
                return var;
            }
        case kwPbState:     return pbProps_ ? pbProps_->pb_state_ : 0;
        case kwPbRMul:      return pbProps_ ? pbProps_->pb_rmul_ : 1.0;
        case kwPbKn:        return pbProps_ ? pbProps_->pb_kn_ : 0;
        case kwPbKs:        return pbProps_ ? pbProps_->pb_ks_ : 0;
		case kwPbKn1:       return pbProps_ ? pbProps_->pb_kn1_ : 0;
		case kwPbKs1:       return pbProps_ ? pbProps_->pb_ks1_ : 0;
        case kwPbMcf:       return pbProps_ ? pbProps_->pb_mcf_ : 1.0;
        case kwPbTStrength: return pbProps_ ? pbProps_->pb_ten_ : 0.0;
		case kwPbCStrength: return pbProps_ ? pbProps_->pb_com_ : 0.0;
		//case kwPbMBStrength: return pbProps_ ? pbProps_->pb_mb_ : 0.0;
		case kwPbTStrength1: return pbProps_ ? pbProps_->pb_ten1_ : 0.0;
		case kwPbCStrength1: return pbProps_ ? pbProps_->pb_com1_ : 0.0;
		//case kwPbMBStrength1: return pbProps_ ? pbProps_->pb_mb1_ : 0.0;
        case kwPbCoh:        return pbProps_ ? pbProps_->pb_coh_ : 0.0;
		case kwPbCoh1:       return pbProps_ ? pbProps_->pb_coh1_ : 0.0;
        //case kwPbFa:         return pbProps_ ? pbProps_->pb_fa_ : 0;
        //case kwPbSig: {
        //       if (!pbProps_ || pbProps_->pb_state_ < 4) return 0.0;
        //        const IContactMechanical *c(convert_getcast<IContactMechanical>(con));
        //        return pbSMax(c).x();
        //    }
        //case kwPbTau: {
        //        if (!pbProps_ || pbProps_->pb_state_ < 4) return 0.0;
        //        const IContactMechanical *c(convert_getcast<IContactMechanical>(con));
        //        return pbSMax(c).y();
        //    }
        case kwPbF: {
                pbProps_ ? var.setValue(pbProps_->pb_F_) : var.setValue(DVect(0.0));
                return var;
            }
        case kwPbM: {
                pbProps_ ? var.setValue(pbProps_->pb_M_) : var.setValue(DAVect(0.0));
                return var;
            }
		case kwMN:			 return pbProps_ ? pbProps_->mn_ : 0.0;
		case kwMS:			 return pbProps_ ? pbProps_->ms_ : 0.0;
		case kwMM:			 return pbProps_ ? pbProps_->mm_ : 0.0;
		case kwFYield:		 return pbProps_ ? pbProps_->Fyield_ : 0.0;
		case kwUCN:		     return pbProps_ ? pbProps_->ucn_ : 0.0;
		case kwUCN1:		 return pbProps_ ? pbProps_->ucn1_ : 0.0;
		case kwUCS:		     return pbProps_ ? pbProps_->ucs_ : 0.0;
		case kwUCB:		     return pbProps_ ? pbProps_->ucb_ : 0.0;
		//case kwALPHA:		 return pbProps_ ? pbProps_->alpha_ : 0.0;
		//case kwTIMING:		 return pbProps_ ? pbProps_->timing_ : 0.0;
		case kwDamagevar:	 return pbProps_ ? pbProps_->D_ : 0.0;
		case kwFy:		     return pbProps_ ? pbProps_->Fy_ : 0.0;
        case kwPbRadius: {
                if (!pbProps_) return 0.0;
                const IContactMechanical *c(convert_getcast<IContactMechanical>(con));
                double Cmax1 = c->getEnd1Curvature().y();
                double Cmax2 = c->getEnd2Curvature().y();
                double br = pbProps_->pb_rmul_ * 1.0 / std::max(Cmax1,Cmax2);
                if (userArea_)
#ifdef THREED
                    br = std::sqrt(userArea_ / dPi);
#else
                    br = userArea_ / 2.0;
#endif
                return br;
            }
        case kwPbEmod: {
                if (!pbProps_) return 0.0;
                const IContactMechanical *c(convert_getcast<IContactMechanical>(con));
                double rsum(0.0);
                if (c->getEnd1Curvature().y())
                    rsum += 1.0/c->getEnd1Curvature().y();
                if (c->getEnd2Curvature().y())
                    rsum += 1.0/c->getEnd2Curvature().y();
                if (userArea_) {
#ifdef THREED
                    double rad = std::sqrt(userArea_ / dPi);
#else
                    double rad = userArea_ / 2.0;
#endif
                    rsum = 2 * rad;
                }
                double emod = pbProps_->pb_kn_ * rsum;
                return emod;
            }
        case kwPbKRatio: {
                if (!pbProps_) return 0.0;
                return (pbProps_->pb_ks_ == 0.0) ? 0.0 : (pbProps_->pb_kn_/pbProps_->pb_ks_);
            }
        }
        assert(0);
        return QVariant();
    }

    bool ContactModelDBondrotation::getPropertyGlobal(uint i) const {
        switch (i) {
        case kwLinF:
        case kwDpF:
        case kwPbF:
            return false;
        }
        return true;
    }

    bool ContactModelDBondrotation::setProperty(uint i,const QVariant &v,IContact *) {
        pbProps pb;
        dpProps dp;

        switch (i) {
        case kwLinKn: {
                if (!v.canConvert<double>())
                    throw Exception("kn must be a double.");
                double val(v.toDouble());
                if (val<0.0)
                    throw Exception("Negative kn not allowed.");
                kn_ = val;
                return true;
            }
        case kwLinKs: {
                if (!v.canConvert<double>())
                    throw Exception("ks must be a double.");
                double val(v.toDouble());
                if (val<0.0)
                    throw Exception("Negative ks not allowed.");
                ks_ = val;
                return true;
            }
        case kwLinFric: {
                if (!v.canConvert<double>())
                    throw Exception("fric must be a double.");
                double val(v.toDouble());
                if (val<0.0)
                    throw Exception("Negative fric not allowed.");
                fric_ = val;
                return false;
            }
        case kwLinF: {
                if (!v.canConvert<DVect>())
                    throw Exception("lin_force must be a vector.");
                DVect val(v.value<DVect>());
                lin_F_ = val;
                return false;
            }
        case kwLinMode: {
                if (!v.canConvert<uint>())
                    throw Exception("lin_mode must be 0 (absolute) or 1 (incremental).");
                uint val(v.toUInt());
                if (val>1)
                    throw Exception("lin_mode must be 0 (absolute) or 1 (incremental).");
                lin_mode_ = val;
                return false;
            }
        case kwRGap: {
                if (!v.canConvert<double>())
                    throw Exception("Reference gap must be a double.");
                double val(v.toDouble());
                rgap_ = val;
                return false;
            }
		case kwNtd: {
			if (!v.canConvert<double>())
				throw Exception("Ntd must be a double.");
			double val(v.toDouble());
			normal_total_disp_ = val;
			return false;
		}
		case kwStdy: {
			if (!v.canConvert<double>())
				throw Exception("Stdy must be a double.");
			double val(v.toDouble());
			shear_total_dispy_ = val;
			return false;
		}
		case kwStdz: {
			if (!v.canConvert<double>())
				throw Exception("Stdz must be a double.");
			double val(v.toDouble());
			shear_total_dispz_ = val;
			return false;
		}
					 
		case kwStd: {
			if (!v.canConvert<double>())
				throw Exception("Std must be a double.");
			double val(v.toDouble());
			shear_total_disp_ = val;
			return false;
		}
		case kwBtry: {
			if (!v.canConvert<double>())
				throw Exception("Btry must be a double.");
			double val(v.toDouble());
			bending_total_roy_ = val;
			return false;
		}
		case kwBtrz: {
			if (!v.canConvert<double>())
				throw Exception("Btrz must be a double.");
			double val(v.toDouble());
			bending_total_roz_ = val;
			return false;
		}
	
		case kwBtr: {
			if (!v.canConvert<double>())
				throw Exception("Btr must be a double.");
			double val(v.toDouble());
			bending_total_ro_ = val;
			return false;
		}
		case kwNpd: {
			if (!v.canConvert<double>())
				throw Exception("Npd must be a double.");
			double val(v.toDouble());
			normal_plastic_disp_ = val;
			return false;
		}
		case kwSpd: {
			if (!v.canConvert<double>())
				throw Exception("Spd must be a double.");
			double val(v.toDouble());
			shear_plastic_disp_ = val;
			return false;
		}
		case kwBpr: {
			if (!v.canConvert<double>())
				throw Exception("Bpr must be a double.");
			double val(v.toDouble());
			bending_plastic_ro_ = val;
			return false;
		}
		/*case kwBendingMy: {
			if (!v.canConvert<double>())
				throw Exception("BendingMy must be a double.");
			double val(v.toDouble());
			BendingMoy_ = val;
			return false;
		}
		case kwBendingMz: {
			if (!v.canConvert<double>())
				throw Exception("BendingMz must be a double.");
			double val(v.toDouble());
			BendingMoz_ = val;
			return false;
		}*/
		/*case kwBendingM: {
			if (!v.canConvert<double>())
				throw Exception("BendingM must be a double.");
			double val(v.toDouble());
			BendingMo_ = val;
			return false;
		}*/
		/*case kwSheary: {
			if (!v.canConvert<double>())
				throw Exception("Sheary must be a double.");
			double val(v.toDouble());
			Vy_ = val;
			return false;
		}
		case kwShearz: {
			if (!v.canConvert<double>())
				throw Exception("Shearz must be a double.");
			double val(v.toDouble());
			Vz_ = val;
			return false;
		}*/




        case kwPbRMul: {
                if (!v.canConvert<double>())
                    throw Exception("pb_rmul must be a double.");
                double val(v.toDouble());
                if (val<=0.0)
                    throw Exception("pb_rmul must be positive.");
                if (val == 1.0 && !pbProps_)
                    return false;
                if (!pbProps_)
                    pbProps_ = NEWC(pbProps());
                pbProps_->pb_rmul_ = val;
                return false;
            }
        case kwPbKn: {
                if (!v.canConvert<double>())
                    throw Exception("pb_kn must be a double.");
                double val(v.toDouble());
                if (val<0.0)
                    throw Exception("Negative pb_kn not allowed.");
                if (val == 0.0 && !pbProps_)
                    return false;
                if (!pbProps_)
                    pbProps_ = NEWC(pbProps());
                pbProps_->pb_kn_ = val;
                return true;
            }
        case kwPbKs: {
                if (!v.canConvert<double>())
                    throw Exception("pb_ks must be a double.");
                double val(v.toDouble());
                if (val<0.0)
                    throw Exception("Negative pb_ks not allowed.");
                if (val == 0.0 && !pbProps_)
                    return false;
                if (!pbProps_)
                    pbProps_ = NEWC(pbProps());
                pbProps_->pb_ks_ = val;
                return true;
            }
		case kwPbKn1: {
			if (!v.canConvert<double>())
				throw Exception("pb_kn1 must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Negative pb_kn1 not allowed.");
			if (val == 0.0 && !pbProps_)
				return false;
			if (!pbProps_)
				pbProps_ = NEWC(pbProps());
			pbProps_->pb_kn1_ = val;
			return true;
		}
		case kwPbKs1: {
			if (!v.canConvert<double>())
				throw Exception("pb_ks1 must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Negative pb_ks1 not allowed.");
			if (val == 0.0 && !pbProps_)
				return false;
			if (!pbProps_)
				pbProps_ = NEWC(pbProps());
			pbProps_->pb_ks1_ = val;
			return true;
		}
        case kwPbMcf: {
                if (!v.canConvert<double>())
                    throw Exception("pb_mcf must be a double.");
                double val(v.toDouble());
                if (val<0.0)
                    throw Exception("Negative pb_mcf not allowed.");
                if (val > 1.0)
                    throw Exception("pb_mcf must be lower or equal to 1.0.");
                if (val == 1.0 && !pbProps_)
                    return false;
                if (!pbProps_)
                    pbProps_ = NEWC(pbProps());
                pbProps_->pb_mcf_ = val;
                return false;
            }
        case kwPbTStrength: {
                if (!v.canConvert<double>())
                    throw Exception("pb_ten must be a double.");
                double val(v.toDouble());
                if (val < 0.0)
                    throw Exception("Negative pb_ten not allowed.");
                if (val == 0.0 && !pbProps_)
                    return false;
                if (!pbProps_)
                    pbProps_ = NEWC(pbProps());
                pbProps_->pb_ten_ = val;
                return false;
            }
		case kwPbCStrength: {
			if (!v.canConvert<double>())
				throw Exception("pb_com must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Negative pb_com not allowed.");
			if (val == 0.0 && !pbProps_)
				return false;
			if (!pbProps_)
				pbProps_ = NEWC(pbProps());
			pbProps_->pb_com_ = val;
			return false;
		}
		/*
		case kwPbMBStrength: {
			if (!v.canConvert<double>())
				throw Exception("pb_mb must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Negative pb_mb not allowed.");
			if (val == 0.0 && !pbProps_)
				return false;
			if (!pbProps_)
				pbProps_ = NEWC(pbProps());
			pbProps_->pb_mb_ = val;
			return false;
		}
		*/
		case kwPbTStrength1: {
			if (!v.canConvert<double>())
				throw Exception("pb_ten1 must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Negative pb_ten1 not allowed.");
			if (val == 0.0 && !pbProps_)
				return false;
			if (!pbProps_)
				pbProps_ = NEWC(pbProps());
			pbProps_->pb_ten1_ = val;
			return false;
		}
		case kwPbCStrength1: {
			if (!v.canConvert<double>())
				throw Exception("pb_com1 must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Negative pb_com1 not allowed.");
			if (val == 0.0 && !pbProps_)
				return false;
			if (!pbProps_)
				pbProps_ = NEWC(pbProps());
			pbProps_->pb_com1_ = val;
			return false;
		}
		/*
		case kwPbMBStrength1: {
			if (!v.canConvert<double>())
				throw Exception("pb_mb1 must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Negative pb_mb1 not allowed.");
			if (val == 0.0 && !pbProps_)
				return false;
			if (!pbProps_)
				pbProps_ = NEWC(pbProps());
			pbProps_->pb_mb1_ = val;
			return false;
		}
	    */
        case kwPbCoh: {
                if (!v.canConvert<double>())
                    throw Exception("pb_coh must be a double.");
                double val(v.toDouble());
                if (val<0.0)
                    throw Exception("Negative pb_coh not allowed.");
                if (val == 0.0 && !pbProps_)
                    return false;
                if (!pbProps_)
                    pbProps_ = NEWC(pbProps());
                pbProps_->pb_coh_ = val;
                return false;
            }
		case kwPbCoh1: {
			if (!v.canConvert<double>())
				throw Exception("pb_coh1 must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Negative pb_coh1 not allowed.");
			if (val == 0.0 && !pbProps_)
				return false;
			if (!pbProps_)
				pbProps_ = NEWC(pbProps());
			pbProps_->pb_coh1_ = val;
			return false;
		}
        //case kwPbFa: {
        //        if (!v.canConvert<double>())
        //            throw Exception("pb_fa must be a double.");
        //        double val(v.toDouble());
        //        if (val<0.0)
        //            throw Exception("Negative pb_fa not allowed.");
        //        if (val>=90.0)
        //            throw Exception("pb_fa must be lower than 90.0 degrees.");
        //        if (val == 0.0 && !pbProps_)
        //            return false;
        //        if (!pbProps_)
        //            pbProps_ = NEWC(pbProps());
        //        pbProps_->pb_fa_ = val;
        //        return false;
        //    }
        case kwPbF: {
                if (!v.canConvert<DVect>())
                    throw Exception("pb_force must be a vector.");
                DVect val(v.value<DVect>());
                if (val.fsum() == 0.0 && !pbProps_)
                    return false;
                if (!pbProps_)
                    pbProps_ = NEWC(pbProps());
                pbProps_->pb_F_ = val;
                return false;
            }
        case kwPbM: {
                DAVect val(0.0);
#ifdef TWOD
                if (!v.canConvert<DAVect>() && !v.canConvert<double>())
                    throw Exception("pb_moment must be an angular vector.");
                if (v.canConvert<DAVect>())
                    val = DAVect(v.value<DAVect>());
                else
                    val = DAVect(v.toDouble());
#else
                if (!v.canConvert<DAVect>() && !v.canConvert<DVect>())
                    throw Exception("pb_moment must be an angular vector.");
                if (v.canConvert<DAVect>())
                    val = DAVect(v.value<DAVect>());
                else
                    val = DAVect(v.value<DVect>());
#endif
                if (val.fsum() == 0.0 && !pbProps_)
                    return false;
                if (!pbProps_)
                    pbProps_ = NEWC(pbProps());
                pbProps_->pb_M_ = val;
                return false;
            }
		case kwMN: {
			if (!v.canConvert<double>())
				throw Exception("mn must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Negative mn not allowed.");
			if (val == 0.0 && !pbProps_)
				return false;
			if (!pbProps_)
				pbProps_ = NEWC(pbProps());
			pbProps_->mn_ = val;
			return false;
		}
		case kwMS: {
			if (!v.canConvert<double>())
				throw Exception("ms must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Negative ms not allowed.");
			if (val == 0.0 && !pbProps_)
				return false;
			if (!pbProps_)
				pbProps_ = NEWC(pbProps());
			pbProps_->ms_ = val;
			return false;
		}
		case kwMM: {
			if (!v.canConvert<double>())
				throw Exception("mm must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Negative mm not allowed.");
			if (val == 0.0 && !pbProps_)
				return false;
			if (!pbProps_)
				pbProps_ = NEWC(pbProps());
			pbProps_->mm_ = val;
			return false;
		}
		case kwFYield: {
			if (!v.canConvert<double>())
				throw Exception("Fyield must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Negative Fyield not allowed.");
			if (val == 0.0 && !pbProps_)
				return false;
			if (!pbProps_)
				pbProps_ = NEWC(pbProps());
			pbProps_->Fyield_ = val;
			return false;
		}
		case kwUCN: {
			if (!v.canConvert<double>())
				throw Exception("ucn must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Negative ucn not allowed.");
			if (val == 0.0 && !pbProps_)
				return false;
			if (!pbProps_)
				pbProps_ = NEWC(pbProps());
			pbProps_->ucn_ = val;
			return false;
		}
		
		case kwUCN1: {
			if (!v.canConvert<double>())
				throw Exception("ucn1 must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Negative ucn1 not allowed.");
			if (val == 0.0 && !pbProps_)
				return false;
			if (!pbProps_)
				pbProps_ = NEWC(pbProps());
			pbProps_->ucn1_ = val;
			return false;
		}
		
		case kwUCS: {
			if (!v.canConvert<double>())
				throw Exception("ucs must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Negative ucs not allowed.");
			if (val == 0.0 && !pbProps_)
				return false;
			if (!pbProps_)
				pbProps_ = NEWC(pbProps());
			pbProps_->ucs_ = val;
			return false;
		}
		case kwUCB: {
			if (!v.canConvert<double>())
				throw Exception("ucb must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Negative ucb not allowed.");
			if (val == 0.0 && !pbProps_)
				return false;
			if (!pbProps_)
				pbProps_ = NEWC(pbProps());
			pbProps_->ucb_ = val;
			return false;
		}


		/*case kwALPHA: {
			if (!v.canConvert<double>())
				throw Exception("alpha must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Negative alpha not allowed.");
			if (val > 1.0)
				throw Exception("alpha must be lower or equal to 1.0.");
			if (val == 1.0 && !pbProps_)
				return false;
			if (!pbProps_)
				pbProps_ = NEWC(pbProps());
			pbProps_->alpha_ = val;
			return false;
		}

		case kwTIMING: {
			if (!v.canConvert<double>())
				throw Exception("timing must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Negative timing not allowed.");
			if (val > 1.0)
				throw Exception("timing must be lower or equal to 1.0.");
			if (val == 1.0 && !pbProps_)
				return false;
			if (!pbProps_)
				pbProps_ = NEWC(pbProps());
			pbProps_->timing_ = val;
			return false;
		}
		*/

		case kwDamagevar: {
			if (!v.canConvert<double>())
				throw Exception("D_ must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Negative D_ not allowed.");
			if (val == 0.0 && !pbProps_)
				return false;
			if (!pbProps_)
				pbProps_ = NEWC(pbProps());
			pbProps_->D_ = val;
			return false;
		}
		case kwFy: {
			if (!v.canConvert<double>())
				throw Exception("Fy_ must be a double.");
			double val(v.toDouble());
			if (val == 0.0 && !pbProps_)
				return false;
			if (!pbProps_)
				pbProps_ = NEWC(pbProps());
			pbProps_->Fy_ = val;
			return false;
		}

        case kwDpNRatio: {
                if (!v.canConvert<double>())
                    throw Exception("dp_nratio must be a double.");
                double val(v.toDouble());
                if (val<0.0)
                    throw Exception("Negative dp_nratio not allowed.");
                if (val == 0.0 && !dpProps_)
                    return false;
                if (!dpProps_)
                    dpProps_ = NEWC(dpProps());
                dpProps_->dp_nratio_ = val;
                return true;
            }
        case kwDpSRatio: {
                if (!v.canConvert<double>())
                    throw Exception("dp_sratio must be a double.");
                double val(v.toDouble());
                if (val<0.0)
                    throw Exception("Negative dp_sratio not allowed.");
                if (val == 0.0 && !dpProps_)
                    return false;
                if (!dpProps_)
                    dpProps_ = NEWC(dpProps());
                dpProps_->dp_sratio_ = val;
                return true;
            }
        case kwDpMode: {
                if (!v.canConvert<int>())
                    throw Exception("The viscous mode dp_mode must be 0, 1, 2, or 3.");
                int val(v.toInt());
                if (val == 0 && !dpProps_)
                    return false;
                if (val < 0 || val > 3)
                    throw Exception("The viscous mode dp_mode must be 0, 1, 2, or 3.");
                if (!dpProps_)
                    dpProps_ = NEWC(dpProps());
                dpProps_->dp_mode_ = val;
                return false;
            }
        case kwDpF: {
                if (!v.canConvert<DVect>())
                    throw Exception("dp_force must be a vector.");
                DVect val(v.value<DVect>());
                if (val.fsum() == 0.0 && !dpProps_)
                    return false;
                if (!dpProps_)
                    dpProps_ = NEWC(dpProps());
                dpProps_->dp_F_ = val;
                return false;
            }
        case kwUserArea: {
                if (!v.canConvert<double>())
                    throw Exception("user_area must be a double.");
                double val(v.toDouble());
                if (val < 0.0)
                    throw Exception("Negative user_area not allowed.");
                userArea_ = val;
                return true;
            }
        }
//    assert(0);
        return false;
    }

    bool ContactModelDBondrotation::getPropertyReadOnly(uint i) const {
        switch (i) {
        case kwDpF:
        case kwLinS:
        case kwEmod:
        case kwKRatio:
        case kwPbState:
        case kwPbRadius:
        //case kwPbSig:
        //case kwPbTau:
        case kwPbEmod:
        case kwPbKRatio:
            return true;
        default:
            break;
        }
        return false;
    }

    bool ContactModelDBondrotation::supportsInheritance(uint i) const {
        switch (i) {
        case kwLinKn:
        case kwLinKs:
        case kwLinFric:
            return true;
        default:
            break;
        }
        return false;
    }

    QString  ContactModelDBondrotation::getMethodArguments(uint i) const {
        QString def = QString();
        switch (i) {
        case kwDeformability:
            return "emod,kratio";
        case kwPbDeformability:
            return "emod,kratio";
        case kwPbBond:
            return "gap";
        case kwPbUnbond:
            return "gap";
        }
        return def;
    }

    bool ContactModelDBondrotation::setMethod(uint i,const QVector<QVariant> &vl,IContact *con) {
        IContactMechanical *c(convert_getcast<IContactMechanical>(con));
        switch (i) {
        case kwDeformability: {
                double emod;
                double krat;
                if (vl.at(0).isNull())
                    throw Exception("Argument emod must be specified with method deformability in contact model %1.",getName());
                emod = vl.at(0).toDouble();
                if (emod<0.0)
                    throw Exception("Negative emod not allowed in contact model %1.",getName());
                if (vl.at(1).isNull())
                    throw Exception("Argument kratio must be specified with method deformability in contact model %1.",getName());
                krat = vl.at(1).toDouble();
                if (krat<0.0)
                    throw Exception("Negative linear stiffness ratio not allowed in contact model %1.",getName());
                double rsq(std::max(c->getEnd1Curvature().y(),c->getEnd2Curvature().y()));
                double rsum(0.0);
                if (c->getEnd1Curvature().y())
                    rsum += 1.0/c->getEnd1Curvature().y();
                if (c->getEnd2Curvature().y())
                    rsum += 1.0/c->getEnd2Curvature().y();
                if (userArea_) {
#ifdef THREED
                    rsq = std::sqrt(userArea_ / dPi);
#else
                    rsq = userArea_ / 2.0;
#endif
                    rsum = rsq + rsq;
                    rsq = 1. / rsq;
                }
#ifdef TWOD
                kn_ = 2.0 * emod / (rsq * rsum);
#else
                kn_ = dPi * emod / (rsq * rsq * rsum);
#endif
                ks_ = (krat == 0.0) ? 0.0 : kn_ / krat;
                setInheritance(1,false);
                setInheritance(2,false);
                return true;
            }
        case kwPbDeformability: {
                //if (!pbProps_ || pbProps_->pb_state_ != 3) return false;
                double emod;
                double krat;
                if (vl.at(0).isNull())
                    throw Exception("Argument emod must be specified with method pb_deformability in contact model %1.",getName());
                emod = vl.at(0).toDouble();
                if (emod<0.0)
                    throw Exception("Negative emod not allowed in contact model %1.",getName());
                if (vl.at(1).isNull())
                    throw Exception("Argument kratio must be specified with method pb_deformability in contact model %1.",getName());
                krat = vl.at(1).toDouble();
                if (krat<0.0)
                    throw Exception("Negative parallel bond stiffness ratio not allowed in contact model %1.",getName());
                double rsum(0.0);
                if (c->getEnd1Curvature().y())
                    rsum += 1.0/c->getEnd1Curvature().y();
                if (c->getEnd2Curvature().y())
                    rsum += 1.0/c->getEnd2Curvature().y();
                if (!pbProps_)
                    pbProps_ = NEWC(pbProps());
                if (userArea_)
#ifdef THREED
                    rsum = 2 * std::sqrt(userArea_ / dPi);
#else
                    rsum = 2 * userArea_ / 2.0;
#endif
                pbProps_->pb_kn_ = emod / rsum;
                pbProps_->pb_ks_ = (krat == 0.0) ? 0.0 : pbProps_->pb_kn_ / krat;
				pbProps_->pb_kn1_ = pbProps_->pb_kn_;         //give pbkn1 and pbks1 an initial value
				pbProps_->pb_ks1_ = pbProps_->pb_ks_;
                return true;
            }
        case kwPbBond: {
                if (pbProps_ && pbProps_->pb_state_ == 4) return false;
                double mingap = -1.0 * limits<double>::max();
                double maxgap = 0;
                if (vl.at(0).canConvert<Double>())
                    maxgap = vl.at(0).toDouble();
                else if (vl.at(0).canConvert<DVect2>()) {
                    DVect2 value = vl.at(0).value<DVect2>();
                    mingap = value.minComp();
                    maxgap = value.maxComp();
                } else if (!vl.at(0).isNull())
                    throw Exception("gap value %1 not recognized in method bond in contact model %2.",vl.at(0),getName());
                double gap = c->getGap();
                if ( gap >= mingap && gap <= maxgap) {
                    if (pbProps_)
                        pbProps_->pb_state_ = 4;
                    else {
                        pbProps_ = NEWC(pbProps());
                        pbProps_->pb_state_ = 4;
                    }
                    return true;
                }
                return false;
            }
        case kwPbUnbond: {
                if (!pbProps_ || pbProps_->pb_state_ == 0) return false;
                double mingap = -1.0 * limits<double>::max();
                double maxgap = 0;
                if (vl.at(0).canConvert<double>())
                    maxgap = vl.at(0).toDouble();
                else if (vl.at(0).canConvert<DVect2>()) {
                    DVect2 value = vl.at(0).value<DVect2>();
                    mingap = value.minComp();
                    maxgap = value.maxComp();
                }
                else if (!vl.at(0).isNull())
                    throw Exception("gap value %1 not recognized in method unbond in contact model %2.",vl.at(0),getName());
                double gap = c->getGap();
                if ( gap >= mingap && gap <= maxgap) {
                    pbProps_->pb_state_ = 0;
                    return true;
                }
                return false;
            }
        case kwArea: {
                if (!userArea_) {
                    double rsq(1./std::max(c->getEnd1Curvature().y(),c->getEnd2Curvature().y()));
#ifdef THREED
                    userArea_ = rsq * rsq * dPi;
#else
                    userArea_ = rsq * 2.0;
#endif
                }
                return true;
            }
        }
        return false;
    }

    double ContactModelDBondrotation::getEnergy(uint i) const {
        double ret(0.0);
        if (!energies_)
            return ret;
        switch (i) {
        case kwEStrain:  return energies_->estrain_;
        case kwESlip:    return energies_->eslip_;
        case kwEDashpot: return energies_->edashpot_;
        case kwEPbStrain:return energies_->epbstrain_;
        }
        assert(0);
        return ret;
    }
	
    bool ContactModelDBondrotation::getEnergyAccumulate(uint i) const {
        switch (i) {
        case kwEStrain:   return false;
        case kwESlip:     return true;
        case kwEDashpot:  return true;
        case kwEPbStrain: return false;
        }
        assert(0);
        return false;
    }

    void ContactModelDBondrotation::setEnergy(uint i,const double &d) {
        if (!energies_) return;
        switch (i) {
        case kwEStrain:   energies_->estrain_  = d; return;
        case kwESlip:     energies_->eslip_    = d; return;
        case kwEDashpot:  energies_->edashpot_ = d; return;
        case kwEPbStrain: energies_->epbstrain_= d; return;
        }
        assert(0);
        return;
    }

    bool ContactModelDBondrotation::validate(ContactModelMechanicalState *state,const double &) {
        assert(state);
        const IContactMechanical *c = state->getMechanicalContact();
        assert(c);

        if (state->trackEnergy_)
            activateEnergy();

        if (inheritanceField_ & linKnMask)
            updateKn(c);
        if (inheritanceField_ & linKsMask)
            updateKs(c);
        if (inheritanceField_ & linFricMask)
            updateFric(c);

        updateEffectiveStiffness(state);
        return checkActivity(state->gap_);
    }

    static const QString knstr("kn");
    bool ContactModelDBondrotation::updateKn(const IContactMechanical *con) {
        assert(con);
        QVariant v1 = con->getEnd1()->getProperty(knstr);
        QVariant v2 = con->getEnd2()->getProperty(knstr);
        if (!v1.isValid() || !v2.isValid())
            return false;
        double kn1 = v1.toDouble();
        double kn2 = v2.toDouble();
        double val = kn_;
        if (kn1 && kn2)
            kn_ = kn1*kn2/(kn1+kn2);
        else if (kn1)
            kn_ = kn1;
        else if (kn2)
            kn_ = kn2;
        return ( (kn_ != val) );
    }

    static const QString ksstr("ks");
    bool ContactModelDBondrotation::updateKs(const IContactMechanical *con) {
        assert(con);
        QVariant v1 = con->getEnd1()->getProperty(ksstr);
        QVariant v2 = con->getEnd2()->getProperty(ksstr);
        if (!v1.isValid() || !v2.isValid())
            return false;
        double ks1 = v1.toDouble();
        double ks2 = v2.toDouble();
        double val = ks_;
        if (ks1 && ks2)
            ks_ = ks1*ks2/(ks1+ks2);
        else if (ks1)
            ks_ = ks1;
        else if (ks2)
            ks_ = ks2;
        return ( (ks_ != val) );
    }

    static const QString fricstr("fric");
    bool ContactModelDBondrotation::updateFric(const IContactMechanical *con) {
        assert(con);
        QVariant v1 = con->getEnd1()->getProperty(fricstr);
        QVariant v2 = con->getEnd2()->getProperty(fricstr);
        if (!v1.isValid() || !v2.isValid())
            return false;
        double fric1 = std::max(0.0,v1.toDouble());
        double fric2 = std::max(0.0,v2.toDouble());
        double val = fric_;
        fric_ = std::min(fric1,fric2);
        return ( (fric_ != val) );
    }

    bool ContactModelDBondrotation::endPropertyUpdated(const QString &name,const IContactMechanical *c) {
        assert(c);
        QStringList availableProperties = getProperties().simplified().replace(" ","").split(",",QString::SkipEmptyParts);
        QRegExp rx(name,Qt::CaseInsensitive);
        int idx = availableProperties.indexOf(rx)+1;
        bool ret=false;

        if (idx<=0)
            return ret;

        switch(idx) {
        case kwLinKn:  { //kn
                if (inheritanceField_ & linKnMask)
                    ret = updateKn(c);
                break;
            }
        case kwLinKs:  { //ks
                if (inheritanceField_ & linKsMask)
                    ret =updateKs(c);
                break;
            }
        case kwLinFric:  { //fric
                if (inheritanceField_ & linFricMask)
                    updateFric(c);
                break;
            }
        }
        return ret;
    }

    void ContactModelDBondrotation::updateEffectiveStiffness(ContactModelMechanicalState *state) {
        DVect2 ret(kn_,ks_);
        // account for viscous damping
        if (dpProps_) {
            DVect2 correct(1.0);
            if (dpProps_->dp_nratio_)
                correct.rx() = sqrt(1.0+dpProps_->dp_nratio_*dpProps_->dp_nratio_) - dpProps_->dp_nratio_;
            if (dpProps_->dp_sratio_)
                correct.ry() = sqrt(1.0+dpProps_->dp_sratio_*dpProps_->dp_sratio_) - dpProps_->dp_sratio_;
            ret /= (correct*correct);
        }
        effectiveTranslationalStiffness_ = ret;
        if (isBonded()) {
            double Cmin1 = state->end1Curvature_.x();
            double Cmax1 = state->end1Curvature_.y();
            double Cmax2 = state->end2Curvature_.y();
            //double dthick = (Cmin1 == 0.0) ? state->end1Extent_.x() : 0.0;
			double dthick = (Cmin1 == 0.0) ? 1.0 : 0.0;
            double br = pbProps_->pb_rmul_ * 1.0 / std::max(Cmax1,Cmax2);
            if (userArea_)
#ifdef THREED
                br = std::sqrt(userArea_ / dPi);
#else
                br = userArea_ / 2.0;
#endif
            double br2 = br*br;
            double pbArea = dthick <= 0.0 ? dPi*br2 : 2.0*br*dthick;
            double bi = dthick <= 0.0 ? 0.25*pbArea*br2 : 2.0*br*br2*dthick/3.0;
            pbProps_->pbTransStiff_.rx() = pbProps_->pb_kn1_*pbArea;
            pbProps_->pbTransStiff_.ry() = pbProps_->pb_ks1_*pbArea;
#if DIM==3
            pbProps_->pbAngStiff_.rx() = pbProps_->pb_ks1_* 2.0 * bi;
            pbProps_->pbAngStiff_.ry() = pbProps_->pb_kn1_* bi;
#endif
            pbProps_->pbAngStiff_.rz() = pbProps_->pb_kn1_* bi;
        }
    }

    double ContactModelDBondrotation::pbStrainEnergy() const {
        double ret(0.0);
        if (pbProps_->pb_kn1_)
            ret = 0.5 * pbProps_->pb_F_.x() * pbProps_->pb_F_.x() / pbProps_->pbTransStiff_.x();
        if (pbProps_->pb_ks1_) {
            DVect tmp = pbProps_->pb_F_;
            tmp.rx() = 0.0;
            double smag2 = tmp.mag2();
            ret += 0.5 * smag2 / pbProps_->pbTransStiff_.y();
        }

#ifdef THREED
        if (pbProps_->pbAngStiff_.x())
            ret += 0.5 * pbProps_->pb_M_.x() * pbProps_->pb_M_.x() / pbProps_->pbAngStiff_.x();
#endif
        if (pbProps_->pbAngStiff_.z()) {
            DAVect tmp = pbProps_->pb_M_;
#ifdef THREED
            tmp.rx() = 0.0;
            double smag2 = tmp.mag2();
#else
            double smag2 = tmp.z() * tmp.z();
#endif
            ret += 0.5 * smag2 / pbProps_->pbAngStiff_.z();
        }
        return ret;
    }

    bool ContactModelDBondrotation::forceDisplacementLaw(ContactModelMechanicalState *state,const double &timestep) {
        assert(state);

        double overlap = rgap_ - state->gap_;
        DVect trans = state->relativeTranslationalIncrement_;
        double correction = 1.0;

        if (state->activated()) {
            if (cmEvents_[fActivated] >= 0) {
				auto c = state->getContact();
				std::vector<fish::Parameter> arg = { fish::Parameter(c->getIThing()) };
				/*FArray<QVariant,2> arg;
                QVariant v;
                IContact * c = const_cast<IContact*>(state->getContact());
                TPtr<IThing> t(c->getIThing());
                v.setValue(t);
                arg.push_back(v);*/
                IFishCallList *fi = const_cast<IFishCallList*>(state->getProgram()->findInterface<IFishCallList>());
                fi->setCMFishCallArguments(c,arg,cmEvents_[fActivated]);
            }
            if (lin_mode_ == 0 && trans.x()) {
                correction = -1.0*overlap / trans.x();
                if (correction < 0)
                    correction = 1.0;
            }
        }

#ifdef THREED
        DVect norm(trans.x(),0.0,0.0);
#else
        DVect norm(trans.x(),0.0);
#endif
        DAVect ang  = state->relativeAngularIncrement_;
		DVect u_s;
		//if (pbProps_->D_ <= 0.6) {
			//DVect u_s = trans;
			//u_s.rx() = 0.0; //these two lines need to be commented if we add lin-force
		//}
		u_s = trans;
		u_s.rx() = 0.0;

		
		DVect ss_F_old = lin_F_;
		// pb_state==0 and pb_state<=4
		if (pbProps_ && pbProps_->pb_state_ == 4) {
			if ((pbProps_->D_ > 0.0) && (overlap > 0.0)) {
				if (lin_mode_ == 0) {
					lin_F_.rx() = overlap * kn_;
				}
				else {
					lin_F_.rx() -= correction * norm.x() * kn_;
				}
				lin_F_.rx() = std::max(0.0, lin_F_.x());

				DVect u_s = trans;
				u_s.rx() = 0.0;
				DVect sforce = lin_F_ - u_s * ks_ * correction;
				sforce.rx() = 0.0;




				// Make sure that the shear force opposses the direction of translation - otherwise you can
				// have strange behavior
				//for (int j=1; j<dim; ++j)
				//{
				//  qDebug()<<sforce.dof(j)<<trans.dof(j);
				//  if (sign<double>(1,sforce.dof(j)) == sign<double>(1,trans.dof(j)))
				//    sforce.rdof(j) = 0.0;
				//}

				// Resolve sliding
				if (state->canFail_) {
					double crit = lin_F_.x() * fric_;
					double sfmag = sforce.mag();
					if (sfmag > crit) {
						double rat = crit / sfmag;
						sforce *= rat;
						if (!lin_S_ && cmEvents_[fSlipChange] >= 0) {
							auto c = state->getContact();
							std::vector<fish::Parameter> arg = { fish::Parameter(c->getIThing()),
																 fish::Parameter() };
							/*FArray<QVariant, 3> arg;
							QVariant p1;
							IContact * c = const_cast<IContact*>(state->getContact());
							TPtr<IThing> t(c->getIThing());
							p1.setValue(t);
							arg.push_back(p1);
							p1.setValue(0);
							arg.push_back(p1);*/
							IFishCallList* fi = const_cast<IFishCallList*>(state->getProgram()->findInterface<IFishCallList>());
							fi->setCMFishCallArguments(c, arg, cmEvents_[fSlipChange]);
						}
						lin_S_ = true;
					}
					else {
						if (lin_S_) {
							if (cmEvents_[fSlipChange] >= 0) {
								auto c = state->getContact();
								std::vector<fish::Parameter> arg = { fish::Parameter(c->getIThing()),
																	 fish::Parameter((qint64)1) };
								/*FArray<QVariant, 3> arg;
								QVariant p1;
								IContact * c = const_cast<IContact*>(state->getContact());
								TPtr<IThing> t(c->getIThing());
								p1.setValue(t);
								arg.push_back(p1);
								p1.setValue(1);
								arg.push_back(p1);*/
								IFishCallList* fi = const_cast<IFishCallList*>(state->getProgram()->findInterface<IFishCallList>());
								fi->setCMFishCallArguments(c, arg, cmEvents_[fSlipChange]);
							}
							lin_S_ = false;
						}
					}
				}

				sforce.rx() = lin_F_.x();
				lin_F_ = sforce;          // total force in linear contact model
			}
		}//no lin_force before bonds break




		if (pbProps_ && pbProps_->pb_state_ < 4) {
			if (overlap > 0.0) {
				if (lin_mode_ == 0) {
					lin_F_.rx() = overlap * kn_;
				}
				else {
					lin_F_.rx() -= correction * norm.x() * kn_;
				}
				lin_F_.rx() = std::max(0.0, lin_F_.x());

				DVect u_s = trans;
				u_s.rx() = 0.0;
				DVect sforce = lin_F_ - u_s * ks_ * correction;
				sforce.rx() = 0.0;




				// Make sure that the shear force opposses the direction of translation - otherwise you can
				// have strange behavior
				//for (int j=1; j<dim; ++j)
				//{
				//  qDebug()<<sforce.dof(j)<<trans.dof(j);
				//  if (sign<double>(1,sforce.dof(j)) == sign<double>(1,trans.dof(j)))
				//    sforce.rdof(j) = 0.0;
				//}

				// Resolve sliding
				if (state->canFail_) {
					double crit = lin_F_.x() * fric_;
					double sfmag = sforce.mag();
					if (sfmag > crit) {
						double rat = crit / sfmag;
						sforce *= rat;
						if (!lin_S_ && cmEvents_[fSlipChange] >= 0) {
							auto c = state->getContact();
							std::vector<fish::Parameter> arg = { fish::Parameter(c->getIThing()),
																 fish::Parameter() };
							/*FArray<QVariant, 3> arg;
							QVariant p1;
							IContact * c = const_cast<IContact*>(state->getContact());
							TPtr<IThing> t(c->getIThing());
							p1.setValue(t);
							arg.push_back(p1);
							p1.setValue(0);
							arg.push_back(p1);*/
							IFishCallList* fi = const_cast<IFishCallList*>(state->getProgram()->findInterface<IFishCallList>());
							fi->setCMFishCallArguments(c, arg, cmEvents_[fSlipChange]);
						}
						lin_S_ = true;
					}
					else {
						if (lin_S_) {
							if (cmEvents_[fSlipChange] >= 0) {
								auto c = state->getContact();
								std::vector<fish::Parameter> arg = { fish::Parameter(c->getIThing()),
																	 fish::Parameter((qint64)1) };
								/*FArray<QVariant, 3> arg;
								QVariant p1;
								IContact * c = const_cast<IContact*>(state->getContact());
								TPtr<IThing> t(c->getIThing());
								p1.setValue(t);
								arg.push_back(p1);
								p1.setValue(1);
								arg.push_back(p1);*/
								IFishCallList* fi = const_cast<IFishCallList*>(state->getProgram()->findInterface<IFishCallList>());
								fi->setCMFishCallArguments(c, arg, cmEvents_[fSlipChange]);
							}
							lin_S_ = false;
						}
					}
				}

				sforce.rx() = lin_F_.x();
				lin_F_ = sforce;          // total force in linear contact model
			}
		}//no lin_force before bonds break








		/*const IContactMechanical* c;
		double rsum(0.0);
		double flag = 0;
		double surface_gap = 0;
		double BL0 = 0;
		
		double r1 = 1.0 / state->end1Curvature_.y();
		double r2 = 1.0 / state->end2Curvature_.y();
		rsum = r1+r2;
		

		if (flag == 0) {
			surface_gap = state->gap_;
			flag = 1;
		}
		
		BL0 = rsum + surface_gap;
		//doube check_dis = alpha * BL; Bond always exist only set its force and momnet as 0
		*/
		

        // Account for dashpot forces
        if (dpProps_) {
            dpProps_->dp_F_.fill(0.0);
            double vcn(0.0), vcs(0.0);
            setDampCoefficients(state->inertialMass_,&vcn,&vcs);
            // Need to change behavior based on the dp_mode
            if (dpProps_->dp_mode_ == 0)  { // Damp all components
                dpProps_->dp_F_ = u_s * (-1.0* vcs) / timestep; // shear component
                dpProps_->dp_F_ -= norm * vcn / timestep;       // normal component
            } else if (dpProps_->dp_mode_ == 1)  { // limit the tensile
                dpProps_->dp_F_ -= norm * vcn / timestep;       // normal component
                if (dpProps_->dp_F_.x() + lin_F_.x() < 0)
                    dpProps_->dp_F_.rx() = - lin_F_.rx();
            } else if (dpProps_->dp_mode_ == 2)  { // limit the shear
                if (!lin_S_)
                    dpProps_->dp_F_ = u_s * (-1.0* vcs) / timestep; // shear component
            } else {
                if (!lin_S_)
                    dpProps_->dp_F_ = u_s * (-1.0* vcs) / timestep; // shear component
                dpProps_->dp_F_ -= norm * vcn / timestep;       // normal component
                if (dpProps_->dp_F_.x() + lin_F_.x() < 0)
                    dpProps_->dp_F_.rx() = - lin_F_.rx();
            }
        }

        // Account for parallel bonds
		//double flag1 = 0;
		//if (flag1==0) {
        bool doPb = false;
        if (pbProps_ && pbProps_->pb_state_ == 4) {
            doPb = true;
            // Parallel Bond Logic:
            // bond area and inertia
            // minimal curvature of end1
            double Cmin1 = state->end1Curvature_.x();
            double Cmax1 = state->end1Curvature_.y();
            double Cmax2 = state->end2Curvature_.y();
            //double dthick = (Cmin1 == 0.0) ? state->end1Extent_.x() : 0.0;
			double dthick = (Cmin1 == 0.0) ? 1.0 : 0.0;
            double br = pbProps_->pb_rmul_ * 1.0 / std::max(Cmax1,Cmax2);
            if (userArea_)
#ifdef THREED
                br = std::sqrt(userArea_ / dPi);
#else
                br = userArea_ / 2.0;
#endif
            double br2 = br*br;
            double pbArea = dthick <= 0.0 ? dPi*br2 : 2.0*br*dthick;
			
            double bi = dthick <= 0.0 ? 0.25*pbArea*br2 : 2.0*br*br2*dthick/3.0;
            pbProps_->pbTransStiff_.rx() = pbProps_->pb_kn1_*pbArea;
            pbProps_->pbTransStiff_.ry() = pbProps_->pb_ks1_*pbArea;

            /* pb force trial */
			double prev_N = pbProps_->pb_F_.x();
            pbProps_->pb_F_ -= norm *(pbProps_->pb_kn1_*pbArea) + u_s * (pbProps_->pb_ks1_*pbArea);

			if (prev_N >= 0.0) {
				if (pbProps_->pb_F_.x() < 0.0) {
					normal_plastic_disp_ = -normal_plastic_disp_;
				}
			}

			if (prev_N < 0.0) {
				if (pbProps_->pb_F_.x() >= 0.0) {
					normal_plastic_disp_ = -normal_plastic_disp_;
				}
			}// 

            /* elastic moment increments */
            //DAVect pbMomentInc(0.0)
#if DIM==3
            pbProps_->pbAngStiff_.rx() = pbProps_->pb_ks1_* 2.0 * bi;
            pbProps_->pbAngStiff_.ry() = pbProps_->pb_kn1_* bi;
#endif
            pbProps_->pbAngStiff_.rz() = pbProps_->pb_kn1_* bi;
			/* pb moment trial */
            DAVect pbMomentInc = ang * pbProps_->pbAngStiff_ *(-1.0);
			DAVect pb_moment_ = pbProps_->pb_M_;
			pbProps_->pb_M_ += pbMomentInc;
			shear_total_dispy_ = shear_total_dispy_+ trans.y();
			shear_total_dispz_ = shear_total_dispz_ + trans.z();
			bending_total_roy_ = bending_total_roy_ + ang.y();
			bending_total_roz_ = bending_total_roz_ + ang.z();
			Vy_ = pbProps_->pb_F_.y();
			Vz_ = pbProps_->pb_F_.z();
			BendingMoy_= pbProps_->pb_M_.y();
			BendingMoz_ = pbProps_->pb_M_.z();

			double BendingMoytilde = pbProps_->pb_M_.y() / br;
			double BendingMoztilde = pbProps_->pb_M_.z() / br;
			double BendingMotilde = std::sqrt(BendingMoztilde*BendingMoztilde + BendingMoytilde * BendingMoytilde);

			DVect bfs(pbProps_->pb_F_);  //force increment
			DAVect bm(pbProps_->pb_M_);  //moment increment
			DAVect bang(ang);            //rotation increment
			bfs.rx() = 0.0;
			bm.rx() = 0.0;
			bang.rx() = 0.0;
			double dbfs = bfs.mag();     //the magnitude of shear force
			double dbm = bm.mag();       //the magnitude of bending moment
			BendingMo_ = dbm;
	        double dnorm = norm.mag();   //the magnitude of the increment of norm disp
			double dtan = u_s.mag();     //the magnitude of the increment of shear disp
			double dbang = bang.mag();   //the magnitude of the increment of bending angle
			double prev_ntd = normal_total_disp_;
			double prev_std = shear_total_disp_;
			double prev_btr = bending_total_ro_;
			shear_total_disp_ = std::sqrt(pow(shear_total_dispy_ ,2.0) + pow(shear_total_dispz_,2.0));
			normal_total_disp_ +=trans.x();
			bending_total_ro_ =std::sqrt(pow(bending_total_roy_, 2.0) + pow(bending_total_roz_, 2.0));

			double Nbar;
			if (pbProps_->pb_F_.x() < 0.0) {
				Nbar = pbProps_->pb_ten1_ * pbArea ;
			}
			else {
			    Nbar = pbProps_->pb_com1_ * pbArea; //compression
			}
			double Vbar = pbProps_->pb_coh1_ * pbArea;
			double Mbar = pbProps_->pb_ten1_ * bi / br / br;

			double Nbar01;
			double Nbar02;
			if (pbProps_->pb_F_.x() < 0.0) {
				Nbar01 = pbProps_->pb_ten_ * pbArea;
			}
			else {
				Nbar02 = pbProps_->pb_com_ * pbArea; //compression
			}
			double Vbar0 = pbProps_->pb_coh_ * pbArea;
			double Mbar0 = pbProps_->pb_ten_ * bi / br / br;


			double N = pbProps_->pb_F_.x();
			double V = dbfs;  //here V and M are the magnitude 
			double M = dbm;

			pbProps_->Fyield_ = pow(std::abs(BendingMotilde / Mbar), 1.001) + pow((N / Nbar), 2.0) + (pow((V / Vbar), 4.0) / (1.0 - pow((N / Nbar), 2.0))) - 1.0;

			//special case because I want the damage can continue
			//if ((pbProps_->D_ > 0.0) && (pbProps_->Fyield_ < 0.0)) {
			//	pbProps_->Fyield_ = 1.0;
			//}
			
			if (pbProps_->Fyield_ > 0.0) {  //find a specific point on the yield surface with shortest dis between this point and the stress point out of the curve
				double NNN;
				double VVV;
				double MMM;
				//the following part is used in one-loading type
				if ((N >= 0.0) && (V==0.0) && (M == 0.0)) {
					NNN = Nbar;
					VVV = 0.0;
				    MMM = 0.0;
				}
				if ((N < 0.0) && (V == 0.0) && (M == 0.0)) {
					NNN = -1.0*Nbar;
					VVV = 0.0;
					MMM = 0.0;
				}
				if ((N == 0.0) && (V != 0.0) && (M == 0.0)) {
					NNN = 0.0;
					VVV = Vbar;
				    MMM = 0.0;
				}
				if ((N == 0.0) && (V == 0.0) && (M != 0.0)) {
					NNN = 0.0;
					VVV = 0.0;
					MMM = Mbar;
				}

				//the following part is used in two-loading type
				double N_trial_new(0.0);
				double V_trial_new(0.0);
				double N_trial_old(0.0);
				double V_trial_old(0.0);
				double M_trial_new(0.0);
				double M_trial_old(0.0);
				double dist_old(1.0e20);
				double dist_new(0.0);

				//for NV
				if ((N != 0.0) && (V != 0.0) && (M == 0.0)) {
					if (N >= 0.0) {
						N_trial_new = (N <= Nbar) ? N : Nbar;
					}
					else {
						N_trial_new = (std::abs(N) <= std::abs(Nbar)) ? N : (-1.0*Nbar);
					}
					//double V_trial_new = 0.0;
					double count_loop = 0.0;
					while (dist_new <= dist_old) {
						N_trial_old = N_trial_new;
						V_trial_old = pow(pow(Vbar,4.0)*pow(1.0 - N_trial_old * N_trial_old / Nbar / Nbar,2.0), 0.25);
						dist_old = pow(N - N_trial_old, 2.0) + pow(V - V_trial_old,2.0);
						if (N >= 0.0) {
							N_trial_new = N_trial_old - Nbar / 1000.0;
						}
						else {
							N_trial_new = N_trial_old + Nbar / 1000.0;
						}
					    V_trial_new = pow(pow(Vbar,4.0)*pow(1.0 - N_trial_new * N_trial_new / Nbar / Nbar, 2.0), 0.25);
						dist_new = pow(N - N_trial_new, 2.0) + pow(V - V_trial_new, 2.0);
						count_loop = count_loop + 1.0;
						if (count_loop == 1000.0) {
							break;
						}
					}

					NNN = N_trial_old;
					VVV = V_trial_old;
					MMM = 0.0;
				}
				//for NM
				if ((N != 0.0) && (V == 0.0) && (M != 0.0)) {
					if (N >= 0.0) {
						N_trial_new = (N <= Nbar) ? N : Nbar;
					}
					else {
						N_trial_new = (std::abs(N) <= std::abs(Nbar)) ? N : (-1.0*Nbar);
					}
					//double M_trial_new = 0.0;
					double count_loop = 0.0;
					while (dist_new <= dist_old) {
						N_trial_old = N_trial_new;
						double pow1 = 1.0 - N_trial_old * N_trial_old / Nbar / Nbar;
						if (pow1 < 0.0) {
							pow1 = 0.0;
						}
						M_trial_old = Mbar * pow(pow1, 1 / 1.001);
						dist_old = pow(N - N_trial_old, 2.0) + pow(BendingMotilde - M_trial_old, 2.0);
						if (N >= 0.0) {
							N_trial_new = N_trial_old - Nbar / 1000.0;
						}
						else {
							N_trial_new = N_trial_old + Nbar /1000.0;
						}
						double pow2 = 1.0 - N_trial_new * N_trial_new / Nbar / Nbar;
						if (pow2 < 0.0) {
							pow2 = 0.0;
						}
						M_trial_new = Mbar * pow(pow2, 1 / 1.001);
						dist_new = pow(N - N_trial_new, 2.0) + pow(BendingMotilde - M_trial_new, 2.0);
						count_loop = count_loop + 1.0;
						if (count_loop == 1000.0) {
							break;
						}
					}

					NNN = N_trial_old;
					VVV = 0.0;
					MMM = M_trial_old;
				}
				//for VM
				if ((N == 0.0) && (V != 0.0) && (M != 0.0)) {
					double V_trial_new = (V <= Vbar) ? V : Vbar;
					double count_loop = 0.0;
					while (dist_new <= dist_old) {
						V_trial_old = V_trial_new;
						double pow3 = 1.0 - V_trial_old * V_trial_old * V_trial_old * V_trial_old / Vbar / Vbar / Vbar / Vbar;
						if (pow3 < 0.0) {
							pow3 = 0.0;
						}
						M_trial_old = Mbar * pow(pow3, 1.0 / 1.001);
						dist_old = pow(V - V_trial_old, 2.0) + pow(BendingMotilde - M_trial_old, 2.0);
						V_trial_new = V_trial_old - Vbar / 1000.0;
						double pow4 = 1.0 - V_trial_new * V_trial_new * V_trial_new * V_trial_new / Vbar / Vbar / Vbar / Vbar;
						if (pow4 < 0.0) {
							pow4 = 0.0;
						}
						M_trial_new = Mbar * pow(pow4, 1.0 / 1.001);
						dist_new = pow(V - V_trial_new, 2.0) + pow(BendingMotilde - M_trial_new, 2.0);
						count_loop = count_loop + 1.0;
						if (count_loop == 1000.0) {
							break;
						}
					}

					NNN = 0.0;
					VVV = V_trial_old;
					MMM = M_trial_old;
				}
				
				//the following part is used for general cases, i.e. NVM

				double temp(0.0);
				double MMM1;
				double MMM2;
				double step1;
				double step2;
				double range1;
				double range2;
				double range3;
				double range4;
				int signal = 0;
				int symbol1 = 0;
				int symbol2 = 0;
				double distance;
				double prev_distance;
				double NL;
				double VL;
				double BendingMotildeL;
				double lambda;
				double prev_NL;
				double prev_VL;
				double prev_BendingMotildeL;
				double prev_lambda;
				
				
				if ((N != 0.0) && (V != 0.0) && (M != 0.0)) {
					if (N >= 0.0) {
						distance = 1.0e200;
						for (int i = 1; i <= 100; i++) {
							if (i == 1) {
								NL = N;
								VL = V;
								BendingMotildeL = BendingMotilde;
								lambda = 0.0;
							}

							else {
								NL = prev_NL;
								VL = prev_VL;
								BendingMotildeL = prev_BendingMotildeL;
								lambda = prev_lambda;
							}
							double g1 = 2.0*(NL - N) + lambda * (2.0*NL / Nbar / Nbar + (2.0*NL / Nbar / Nbar * pow(VL, 4.0) / pow(Vbar, 4.0)) / (1.0 - NL * NL / Nbar / Nbar) / (1.0 - NL * NL / Nbar / Nbar));
							double g2 = 2.0*(VL - V) + lambda * (4 * VL*VL*VL / (pow(Vbar, 4.0)) / (1.0 - NL * NL / Nbar / Nbar));
							double g3 = 2.0*(BendingMotildeL - BendingMotilde) + lambda * (1.001*pow(BendingMotildeL / Mbar, 0.001) / Mbar);
							double phi = pow(BendingMotildeL / Mbar, 1.001) + pow(NL / Nbar, 2.0) + pow(VL / Vbar, 4.0) / (1.0 - pow(NL / Nbar, 2.0)) - 1.0;

							double Mat11 = 2.0 + lambda * (2.0 / Nbar / Nbar + (2.0 * pow(VL, 4.0) *pow(1.0 - NL * NL / Nbar / Nbar, 2.0) / Nbar / Nbar / (pow(Vbar, 4.0)) + 8.0*NL*NL*pow(VL, 4.0) / (pow(Nbar, 4.0)) / (pow(Vbar, 4.0))* (1.0 - NL * NL / Nbar / Nbar)) / (pow(1.0 - NL * NL / Nbar / Nbar, 4.0)));
							double Mat12 = lambda * (8.0 * NL*VL*VL*VL / (pow(Vbar, 4.0)) / Nbar / Nbar) / (pow(1.0 - NL * NL / Nbar / Nbar, 2.0));
							double Mat13 = 0.0;
							double Mat14 = 2.0*NL / Nbar / Nbar + (2 * NL*pow(VL, 4.0) / (pow(Vbar, 4.0)) / Nbar / Nbar) / (pow(1.0 - NL * NL / Nbar / Nbar, 2.0));

							double Mat21 = lambda * (8.0 * NL*pow(VL, 3.0) / (pow(Vbar, 4.0)) / Nbar / Nbar) / (pow(1.0 - NL * NL / Nbar / Nbar, 2.0));
							double Mat22 = 2.0 + lambda * (12.0 * VL*VL / (pow(Vbar, 4.0))) / (1.0 - NL * NL / Nbar / Nbar);
							double Mat23 = 0.0;
							double Mat24 = 4.0*VL*VL*VL / (pow(Vbar, 4.0)) / (1.0 - NL * NL / Nbar / Nbar);

							double Mat31 = 0.0;
							double Mat32 = 0.0;
							double Mat33 = 2.0 + lambda * (1.001*0.001 / BendingMotildeL / BendingMotildeL * pow(BendingMotildeL / Mbar, -0.999));
							double Mat34 = pow(BendingMotildeL / Mbar, 0.001) * 1.001 / Mbar;

							double Mat41 = 2.0*NL / Nbar / Nbar + (2.0 * NL*pow(VL, 4.0) / (pow(Vbar, 4.0) / Nbar / Nbar)) / (pow(1.0 - NL * NL / Nbar / Nbar, 2.0));
							double Mat42 = 4.0*VL*VL*VL / (pow(Vbar, 4.0)) / (1.0 - NL * NL / Nbar / Nbar);
							double Mat43 = 1.001*pow(BendingMotildeL / Mbar, 0.001) / Mbar;
							double Mat44 = 0.0;
							double Mat[4][4] = { {Mat11,Mat12,Mat13,Mat14},{Mat21,Mat22,Mat23,Mat24},{Mat31,Mat32,Mat33,Mat34},{Mat41,Mat42,Mat43,Mat44} };

							double m0 = Mat11; double m1 = Mat12; double m2 = Mat13; double m3 = Mat14;
							double m4 = Mat21; double m5 = Mat22; double m6 = Mat23; double m7 = Mat24;
							double m8 = Mat31; double m9 = Mat32; double m10 = Mat33; double m11 = Mat34;
							double m12 = Mat41; double m13 = Mat42; double m14 = Mat43; double m15 = Mat44;

							double inv0 = m5 * m10 * m15 - m5 * m11  * m14 - m9 * m6   * m15 + m9 * m7   * m14 + m13 * m6   * m11 - m13 * m7   * m10;
							double inv4 = -m4 * m10  * m15 + m4 * m11  * m14 + m8 * m6   * m15 - m8 * m7   * m14 - m12 * m6   * m11 + m12 * m7   * m10;
							double inv8 = m4 * m9  * m15 - m4 * m11  * m13 - m8 * m5  * m15 + m8 * m7  * m13 + m12 * m5  * m11 - m12 * m7  * m9;
							double inv12 = -m4 * m9  * m14 + m4 * m10  * m13 + m8 * m5  * m14 - m8 * m6  * m13 - m12 * m5  * m10 + m12 * m6  * m9;

							double inv1 = -m1 * m10 * m15 + m1 * m11  * m14 + m9 * m2  * m15 - m9 * m3  * m14 - m13 * m2  * m11 + m13 * m3  * m10;
							double inv5 = m0 * m10  * m15 - m0 * m11  * m14 - m8 * m2  * m15 + m8 * m3  * m14 + m12 * m2  * m11 - m12 * m3  * m10;
							double inv9 = -m0 * m9  * m15 + m0 * m11  * m13 + m8 * m1  * m15 - m8 * m3  * m13 - m12 * m1  * m11 + m12 * m3  * m9;
							double inv13 = m0 * m9  * m14 - m0 * m10  * m13 - m8 * m1  * m14 + m8 * m2  * m13 + m12 * m1  * m10 - m12 * m2  * m9;

							double inv2 = m1 * m6  * m15 - m1 * m7  * m14 - m5 * m2  * m15 + m5 * m3  * m14 + m13 * m2  * m7 - m13 * m3  * m6;
							double inv6 = -m0 * m6  * m15 + m0 * m7  * m14 + m4 * m2  * m15 - m4 * m3  * m14 - m12 * m2  * m7 + m12 * m3  * m6;
							double inv10 = m0 * m5  * m15 - m0 * m7  * m13 - m4 * m1  * m15 + m4 * m3  * m13 + m12 * m1  * m7 - m12 * m3  * m5;
							double inv14 = -m0 * m5  * m14 + m0 * m6  * m13 + m4 * m1  * m14 - m4 * m2  * m13 - m12 * m1  * m6 + m12 * m2  * m5;

							double inv3 = -m1 * m6  * m11 + m1 * m7  * m10 + m5 * m2  * m11 - m5 * m3  * m10 - m9 * m2  * m7 + m9 * m3  * m6;
							double inv7 = m0 * m6  * m11 - m0 * m7  * m10 - m4 * m2  * m11 + m4 * m3  * m10 + m8 * m2  * m7 - m8 * m3  * m6;
							double inv11 = -m0 * m5  * m11 + m0 * m7  * m9 + m4 * m1  * m11 - m4 * m3  * m9 - m8 * m1  * m7 + m8 * m3  * m5;
							double inv15 = m0 * m5  * m10 - m0 * m6  * m9 - m4 * m1  * m10 + m4 * m2  * m9 + m8 * m1  * m6 - m8 * m2  * m5;

							double det = m0 * inv0 + m1 * inv4 + m2 * inv8 + m3 * inv12;
							if (det > 0.0) {
								det = 1.0 / det;
							}
							double inverse_Mat11 = inv0 * det; double inverse_Mat12 = inv1 * det; double inverse_Mat13 = inv2 * det; double inverse_Mat14 = inv3 * det;
							double inverse_Mat21 = inv4 * det; double inverse_Mat22 = inv5 * det; double inverse_Mat23 = inv6 * det; double inverse_Mat24 = inv7 * det;
							double inverse_Mat31 = inv8 * det; double inverse_Mat32 = inv9 * det; double inverse_Mat33 = inv10 * det; double inverse_Mat34 = inv11 * det;
							double inverse_Mat41 = inv12 * det; double inverse_Mat42 = inv13 * det; double inverse_Mat43 = inv14 * det; double inverse_Mat44 = inv15 * det;
							prev_NL = NL; prev_VL = VL; prev_BendingMotildeL = BendingMotildeL; prev_lambda = lambda;
							NL = NL - (inverse_Mat11*g1 + inverse_Mat12 * g2 + inverse_Mat13 * g3 + inverse_Mat14 * phi);
							VL = VL - (inverse_Mat21*g1 + inverse_Mat22 * g2 + inverse_Mat23 * g3 + inverse_Mat24 * phi);
							BendingMotildeL = BendingMotildeL - (inverse_Mat31*g1 + inverse_Mat32 * g2 + inverse_Mat33 * g3 + inverse_Mat34 * phi);

							if (i > 1) {
								if (NL < 0.0 || VL < 0.0 || BendingMotildeL < 0.0) {
									NNN = prev_NL;
									VVV = prev_VL;
									MMM = prev_BendingMotildeL;
									break;
								}
							}
							if (i==1) {
								if (NL < 0.0 || VL < 0.0 || BendingMotildeL < 0.0) {
									signal = 1;
									break;
								}
							}

							lambda = lambda - (inverse_Mat41*g1 + inverse_Mat42 * g2 + inverse_Mat43 * g3 + inverse_Mat44 * phi);
							prev_distance = distance;
							distance = pow(NL - N, 2.0) + pow(VL - V, 2.0) + pow(BendingMotildeL - BendingMotilde, 2.0);

							if (i > 1) {
								if (distance > prev_distance) {
									NNN = prev_NL;
									VVV = prev_VL;
									MMM = prev_BendingMotildeL;
									break;
								}
							}

							if (distance < 1.0e-3) {
								NNN = prev_NL;
								VVV = prev_VL;
								MMM = prev_BendingMotildeL;
								break;
							}

							if (i == 100) {
								NNN = NL;
								VVV = VL;
								MMM = BendingMotildeL;
								break;
							}
						}//newton method
						
						if (signal == 1) {// some special cases due to the form of the yield surface
							double N_search_up = (N <= Nbar) ? N : (0.999999*Nbar);
							double N_search_low = 0.8 * N_search_up;
							double V_search_up = (V <= Vbar) ? V : Vbar;
							double V_search_low = 0.5 * V_search_up;
							range1 = (2.0 / 3.0)*Nbar02;
							range2 = (1.0 / 3.0)*Nbar02;
							range3 = (2.0 / 3.0)*Vbar0;
							range4 = (1.0 / 3.0)*Vbar0;
							if (range1 < N) {
								step1 = (N_search_up - N_search_low) / 99.0;
							}
							if (range3 < V) {
								step2 = (V_search_up - V_search_low) / 49.0;
							}
							if ((range2 < N) & (N <= range1)) {
								step1 = (N_search_up - N_search_low) / 69.0;
							}
							if ((range4 < V) & (V <= range3)) {
								step2 = (V_search_up - V_search_low) / 39.0;
							}
							if (N <= range2) {
								step1 = (N_search_up - N_search_low) / 49.0;
							}
							if (V <= range4) {
								step2 = (V_search_up - V_search_low) / 29.0;
							}



							double Nini;
							double Vini;
							double distance;
							double min_distance = 1.0e200;
							for (Nini = N_search_low; Nini <= N_search_up; Nini = Nini + step1) {
								for (Vini = V_search_low; Vini <= V_search_up; Vini = Vini + step2) {
									temp = 1.0 - pow(Nini / Nbar, 2.0) - pow(Vini / Vbar, 4.0) / (1.0 - pow(Nini / Nbar, 2.0));
									if (temp >= 0.0) {
										symbol1 = 1;
										MMM1 = Mbar * pow(temp, 1.0 / 1.001);
										if (MMM1 < Mbar) {
											symbol2 = 1;
											MMM2 = MMM1;
											distance = pow(N - Nini, 2.0) + pow(V - Vini, 2.0) + pow(BendingMotilde - MMM2, 2.0);
											if (distance < min_distance) {
												min_distance = distance;
												NNN = Nini;
												VVV = Vini;
												MMM = MMM2;
												if (min_distance < 1e-2) {
													break;
													NNN = Nini;
													VVV = Vini;
													MMM = MMM2;
												}
											}
										}
									}
								}
							}
							if (symbol1 == 0 || symbol2 == 0) {
								NNN = 0.0;
								VVV = 0.0;
								MMM = 0.0;
								pbProps_->D_ = 0.999;;
							}

						}

					}//N>0
					else {
						distance = 1.0e200;
						for (int i = 1; i <= 100; i++) {
							if (i == 1) {
								NL = N;
								VL = V;
								BendingMotildeL = BendingMotilde;
								lambda = 0.0;
							}

							else {
								NL = prev_NL;
								VL = prev_VL;
								BendingMotildeL = prev_BendingMotildeL;
								lambda = prev_lambda;
							}
							double g1 = 2.0*(NL - N) + lambda * (2.0*NL / Nbar / Nbar + (2.0*NL / Nbar / Nbar * pow(VL, 4.0) / pow(Vbar, 4.0)) / (1.0 - NL * NL / Nbar / Nbar) / (1.0 - NL * NL / Nbar / Nbar));
							double g2 = 2.0*(VL - V) + lambda * (4 * VL*VL*VL / (pow(Vbar, 4.0)) / (1.0 - NL * NL / Nbar / Nbar));
							double g3 = 2.0*(BendingMotildeL - BendingMotilde) + lambda * (1.001*pow(BendingMotildeL / Mbar, 0.001) / Mbar);
							double phi = pow(BendingMotildeL / Mbar, 1.001) + pow(NL / Nbar, 2.0) + pow(VL / Vbar, 4.0) / (1.0 - pow(NL / Nbar, 2.0)) - 1.0;

							double Mat11 = 2.0 + lambda * (2.0 / Nbar / Nbar + (2.0 * pow(VL, 4.0) *pow(1.0 - NL * NL / Nbar / Nbar, 2.0) / Nbar / Nbar / (pow(Vbar, 4.0)) + 8.0*NL*NL*pow(VL, 4.0) / (pow(Nbar, 4.0)) / (pow(Vbar, 4.0))* (1.0 - NL * NL / Nbar / Nbar)) / (pow(1.0 - NL * NL / Nbar / Nbar, 4.0)));
							double Mat12 = lambda * (8.0 * NL*VL*VL*VL / (pow(Vbar, 4.0)) / Nbar / Nbar) / (pow(1.0 - NL * NL / Nbar / Nbar, 2.0));
							double Mat13 = 0.0;
							double Mat14 = 2.0*NL / Nbar / Nbar + (2 * NL*pow(VL, 4.0) / (pow(Vbar, 4.0)) / Nbar / Nbar) / (pow(1.0 - NL * NL / Nbar / Nbar, 2.0));

							double Mat21 = lambda * (8.0 * NL*pow(VL, 3.0) / (pow(Vbar, 4.0)) / Nbar / Nbar) / (pow(1.0 - NL * NL / Nbar / Nbar, 2.0));
							double Mat22 = 2.0 + lambda * (12.0 * VL*VL / (pow(Vbar, 4.0))) / (1.0 - NL * NL / Nbar / Nbar);
							double Mat23 = 0.0;
							double Mat24 = 4.0*VL*VL*VL / (pow(Vbar, 4.0)) / (1.0 - NL * NL / Nbar / Nbar);

							double Mat31 = 0.0;
							double Mat32 = 0.0;
							double Mat33 = 2.0 + lambda * (1.001*0.001 / BendingMotildeL / BendingMotildeL * pow(BendingMotildeL / Mbar, -0.999));
							double Mat34 = pow(BendingMotildeL / Mbar, 0.001) * 1.001 / Mbar;

							double Mat41 = 2.0*NL / Nbar / Nbar + (2.0 * NL*pow(VL, 4.0) / (pow(Vbar, 4.0) / Nbar / Nbar)) / (pow(1.0 - NL * NL / Nbar / Nbar, 2.0));
							double Mat42 = 4.0*VL*VL*VL / (pow(Vbar, 4.0)) / (1.0 - NL * NL / Nbar / Nbar);
							double Mat43 = 1.001*pow(BendingMotildeL / Mbar, 0.001) / Mbar;
							double Mat44 = 0.0;
							double Mat[4][4] = { {Mat11,Mat12,Mat13,Mat14},{Mat21,Mat22,Mat23,Mat24},{Mat31,Mat32,Mat33,Mat34},{Mat41,Mat42,Mat43,Mat44} };

							double m0 = Mat11; double m1 = Mat12; double m2 = Mat13; double m3 = Mat14;
							double m4 = Mat21; double m5 = Mat22; double m6 = Mat23; double m7 = Mat24;
							double m8 = Mat31; double m9 = Mat32; double m10 = Mat33; double m11 = Mat34;
							double m12 = Mat41; double m13 = Mat42; double m14 = Mat43; double m15 = Mat44;

							double inv0 = m5 * m10 * m15 - m5 * m11  * m14 - m9 * m6   * m15 + m9 * m7   * m14 + m13 * m6   * m11 - m13 * m7   * m10;
							double inv4 = -m4 * m10  * m15 + m4 * m11  * m14 + m8 * m6   * m15 - m8 * m7   * m14 - m12 * m6   * m11 + m12 * m7   * m10;
							double inv8 = m4 * m9  * m15 - m4 * m11  * m13 - m8 * m5  * m15 + m8 * m7  * m13 + m12 * m5  * m11 - m12 * m7  * m9;
							double inv12 = -m4 * m9  * m14 + m4 * m10  * m13 + m8 * m5  * m14 - m8 * m6  * m13 - m12 * m5  * m10 + m12 * m6  * m9;

							double inv1 = -m1 * m10 * m15 + m1 * m11  * m14 + m9 * m2  * m15 - m9 * m3  * m14 - m13 * m2  * m11 + m13 * m3  * m10;
							double inv5 = m0 * m10  * m15 - m0 * m11  * m14 - m8 * m2  * m15 + m8 * m3  * m14 + m12 * m2  * m11 - m12 * m3  * m10;
							double inv9 = -m0 * m9  * m15 + m0 * m11  * m13 + m8 * m1  * m15 - m8 * m3  * m13 - m12 * m1  * m11 + m12 * m3  * m9;
							double inv13 = m0 * m9  * m14 - m0 * m10  * m13 - m8 * m1  * m14 + m8 * m2  * m13 + m12 * m1  * m10 - m12 * m2  * m9;

							double inv2 = m1 * m6  * m15 - m1 * m7  * m14 - m5 * m2  * m15 + m5 * m3  * m14 + m13 * m2  * m7 - m13 * m3  * m6;
							double inv6 = -m0 * m6  * m15 + m0 * m7  * m14 + m4 * m2  * m15 - m4 * m3  * m14 - m12 * m2  * m7 + m12 * m3  * m6;
							double inv10 = m0 * m5  * m15 - m0 * m7  * m13 - m4 * m1  * m15 + m4 * m3  * m13 + m12 * m1  * m7 - m12 * m3  * m5;
							double inv14 = -m0 * m5  * m14 + m0 * m6  * m13 + m4 * m1  * m14 - m4 * m2  * m13 - m12 * m1  * m6 + m12 * m2  * m5;

							double inv3 = -m1 * m6  * m11 + m1 * m7  * m10 + m5 * m2  * m11 - m5 * m3  * m10 - m9 * m2  * m7 + m9 * m3  * m6;
							double inv7 = m0 * m6  * m11 - m0 * m7  * m10 - m4 * m2  * m11 + m4 * m3  * m10 + m8 * m2  * m7 - m8 * m3  * m6;
							double inv11 = -m0 * m5  * m11 + m0 * m7  * m9 + m4 * m1  * m11 - m4 * m3  * m9 - m8 * m1  * m7 + m8 * m3  * m5;
							double inv15 = m0 * m5  * m10 - m0 * m6  * m9 - m4 * m1  * m10 + m4 * m2  * m9 + m8 * m1  * m6 - m8 * m2  * m5;

							double det = m0 * inv0 + m1 * inv4 + m2 * inv8 + m3 * inv12;
							if (det > 0.0) {
								det = 1.0 / det;
							}
							double inverse_Mat11 = inv0 * det; double inverse_Mat12 = inv1 * det; double inverse_Mat13 = inv2 * det; double inverse_Mat14 = inv3 * det;
							double inverse_Mat21 = inv4 * det; double inverse_Mat22 = inv5 * det; double inverse_Mat23 = inv6 * det; double inverse_Mat24 = inv7 * det;
							double inverse_Mat31 = inv8 * det; double inverse_Mat32 = inv9 * det; double inverse_Mat33 = inv10 * det; double inverse_Mat34 = inv11 * det;
							double inverse_Mat41 = inv12 * det; double inverse_Mat42 = inv13 * det; double inverse_Mat43 = inv14 * det; double inverse_Mat44 = inv15 * det;
							prev_NL = NL; prev_VL = VL; prev_BendingMotildeL = BendingMotildeL; prev_lambda = lambda;
							NL = NL - (inverse_Mat11*g1 + inverse_Mat12 * g2 + inverse_Mat13 * g3 + inverse_Mat14 * phi);
							VL = VL - (inverse_Mat21*g1 + inverse_Mat22 * g2 + inverse_Mat23 * g3 + inverse_Mat24 * phi);
							BendingMotildeL = BendingMotildeL - (inverse_Mat31*g1 + inverse_Mat32 * g2 + inverse_Mat33 * g3 + inverse_Mat34 * phi);
							if (i > 1) {
								if (NL > 0.0 || VL < 0.0 || BendingMotildeL < 0.0) {
									NNN = prev_NL;
									VVV = prev_VL;
									MMM = prev_BendingMotildeL;
									break;
								}
							}
							if (i==1) {
								if (NL > 0.0 || VL < 0.0 || BendingMotildeL < 0.0) {
									signal = 1;
									break;
								}
							}

							lambda = lambda - (inverse_Mat41*g1 + inverse_Mat42 * g2 + inverse_Mat43 * g3 + inverse_Mat44 * phi);
							prev_distance = distance;
							distance = pow(NL - N, 2.0) + pow(VL - V, 2.0) + pow(BendingMotildeL - BendingMotilde, 2.0);

							if (i > 1) {
								if (distance > prev_distance) {
									NNN = prev_NL;
									VVV = prev_VL;
									MMM = prev_BendingMotildeL;
									break;
								}
							}

							if (distance < 1.0e-3) {
								NNN = prev_NL;
								VVV = prev_VL;
								MMM = prev_BendingMotildeL;
								break;
							}

							if (i == 100) {
								NNN = NL;
								VVV = VL;
								MMM = BendingMotildeL;
								break;
							}
						}//newton method
						if (signal == 1) {
							double N_search_low = (std::abs(N) <= std::abs(Nbar)) ? N : (-1.0*0.999999*Nbar);
							double N_search_up = 0.5 * N_search_low;
							double V_search_up = (V <= Vbar) ? V : Vbar;
							double V_search_low = 0.5 * V_search_up;
							range1 = (2.0 / 3.0)*Nbar01;
							range2 = (1.0 / 3.0)*Nbar01;
							range3 = (2.0 / 3.0)*Vbar0;
							range4 = (1.0 / 3.0)*Vbar0;
							if (std::abs(N) > range1) {
								step1 = (N_search_up - N_search_low) / 49.0; //positive
							}
							if (range3 < V) {
								step2 = (V_search_up - V_search_low) / 49.0;
							}
							if ((range1 >= std::abs(N)) & (std::abs(N) > range2)) {
								step1 = (N_search_up - N_search_low) / 39.0; //positive
							}
							if ((range4 < V) & (V <= range3)) {
								step2 = (V_search_up - V_search_low) / 39.0;
							}
							if (range2 >= std::abs(N)) {
								step1 = (N_search_up - N_search_low) / 29.0; //positive
							}
							if (V <= range4) {
								step2 = (V_search_up - V_search_low) / 29.0;
							}


							double Nini;
							double Vini;
							double distance;
							double min_distance = 1.0e200;
							for (Nini = N_search_low; Nini <= N_search_up; Nini = Nini + step1) {
								for (Vini = V_search_low; Vini <= V_search_up; Vini = Vini + step2) {
									temp = 1.0 - pow(Nini / Nbar, 2.0) - pow(Vini / Vbar, 4.0) / (1.0 - pow(Nini / Nbar, 2.0));
									if (temp >= 0.0) {
										symbol1 = 1;
										MMM1 = Mbar * pow(temp, 1.0 / 1.001);
										if (MMM1 < Mbar) {
											symbol2 = 1;
											MMM2 = MMM1;
											distance = pow(N - Nini, 2.0) + pow(V - Vini, 2.0) + pow(BendingMotilde - MMM2, 2.0);
											if (distance < min_distance) {
												min_distance = distance;
												NNN = Nini;
												VVV = Vini;
												MMM = MMM2;
												if (min_distance < 1e-2) {
													break;
													NNN = Nini;
													VVV = Vini;
													MMM = MMM2;
												}
											}
										}
									}
								}
							}
							if (symbol1 == 0 || symbol2 == 0) {
								NNN = 0.0;
								VVV = 0.0;
								MMM = 0.0;
								pbProps_->D_ = 0.999;
							}

						}
					}

				}//for general cases
				//here the closest point is found, then the plastic multiplier can be calculated. 
				//the following part is used to calculate the plastic multiplier
				double df_dN = 2.0 * N / Nbar / Nbar + (2.0 * N * pow(V, 4.0) / Nbar / Nbar / pow(Vbar, 4.0)) / (1.0 - N * N / Nbar / Nbar) / (1.0 - N * N / Nbar / Nbar);
				double df_dV = (4.0 * pow(V, 3.0) / pow(Vbar, 4.0)) / (1.0 - N * N / Nbar / Nbar);
				double sign_M = (M > 0.0) - (M < 0.0);
				if (M == 0.0) {
					sign_M = 0.0;
				} 
				double df_dM = 1.001 * (pow(std::abs(M) / Mbar, 0.001))*sign_M*(1.0 / Mbar);

				double dN_dunp = -1.0 * (1.0 - pbProps_->D_) * pbProps_->pb_kn_ * pbArea;
				double dV_dusp = -1.0 * (1.0 - pbProps_->D_) * pbProps_->pb_ks_ * pbArea;
				double dM_dthetabp = -1.0 * (1.0 - pbProps_->D_) * pbProps_->pb_kn_ * bi / br;
				
				double dN_dD = 1.0* pbProps_->pb_kn_ * (normal_total_disp_ - normal_plastic_disp_) * pbArea;
				double dV_dD = -1.0 * pbProps_->pb_ks_ * (shear_total_disp_ - shear_plastic_disp_) * pbArea;
				double dM_dD = -1.0 * pbProps_->pb_kn_ * bi * (bending_total_ro_ - bending_plastic_ro_) / br;
				double sign_Npd = (normal_plastic_disp_ > 0.0) - (normal_plastic_disp_ < 0.0);
				double sign_Spd = (shear_plastic_disp_ > 0.0) - (shear_plastic_disp_ < 0.0);
				double sign_Bpr = (bending_plastic_ro_ > 0.0) - (bending_plastic_ro_ <  0.0);
				if (normal_plastic_disp_ == 0.0) {
					double sign_Npd = 0.0;
				}
				if (shear_plastic_disp_ == 0.0) {
					double sign_Spd = 0.0;
				}
				if (bending_plastic_ro_ == 0.0) {
					double sign_Bpr = 0.0;
				}
				double dD_dunp;
				double dD_dusp;
				double dD_dthetabp;
				if (pbProps_->pb_F_.x() < 0.0) { //tension
					dD_dunp = pow(2.718281828459, (-(std::abs(normal_plastic_disp_) / pbProps_->ucn_ + std::abs(shear_plastic_disp_) / pbProps_->ucs_ + std::abs(bending_plastic_ro_) / pbProps_->ucb_))) / pbProps_->ucn_ * sign_Npd;
					dD_dusp = pow(2.718281828459, (-(std::abs(normal_plastic_disp_) / pbProps_->ucn_ + std::abs(shear_plastic_disp_) / pbProps_->ucs_ + std::abs(bending_plastic_ro_) / pbProps_->ucb_))) / pbProps_->ucs_ * sign_Spd;
					dD_dthetabp = pow(2.718281828459, (-(std::abs(normal_plastic_disp_) / pbProps_->ucn_ + std::abs(shear_plastic_disp_) / pbProps_->ucs_ + std::abs(bending_plastic_ro_) / pbProps_->ucb_))) / pbProps_->ucb_ * sign_Bpr;
				}
				
				else {
					dD_dunp = pow(2.718281828459, (-(std::abs(normal_plastic_disp_) / pbProps_->ucn1_ + std::abs(shear_plastic_disp_) / pbProps_->ucs_ + std::abs(bending_plastic_ro_) / pbProps_->ucb_))) / pbProps_->ucn1_ * sign_Npd;
					dD_dusp = pow(2.718281828459, (-(std::abs(normal_plastic_disp_) / pbProps_->ucn1_ + std::abs(shear_plastic_disp_) / pbProps_->ucs_ + std::abs(bending_plastic_ro_) / pbProps_->ucb_))) / pbProps_->ucs_ * sign_Spd;
					dD_dthetabp = pow(2.718281828459, (-(std::abs(normal_plastic_disp_) / pbProps_->ucn1_ + std::abs(shear_plastic_disp_) / pbProps_->ucs_ + std::abs(bending_plastic_ro_) / pbProps_->ucb_))) / pbProps_->ucb_ * sign_Bpr;
				}
				

				if (std::abs(NNN) == std::abs(Nbar)) {
					if (N >= 0.0) {
						NNN = 0.999999*Nbar;
					}
					else {
						NNN = -0.999999*Nbar;
					}
				}

				pbProps_->mn_ = 2.0 * NNN / Nbar / Nbar + (2.0 * NNN * pow(VVV, 4.0) / Nbar / Nbar / pow(Vbar, 4.0)) / (1.0 - NNN * NNN / Nbar / Nbar) / (1.0 - NNN * NNN / Nbar / Nbar);
				pbProps_->ms_ = (4.0 * pow(VVV, 3.0) / pow(Vbar, 4.0)) / (1.0 - NNN * NNN / Nbar / Nbar);
				double sign_MMM = (MMM > 0.0) - (MMM < 0.0);
				if (MMM == 0.0) {
					sign_MMM = 0.0;
				}
				pbProps_->mm_ = 1.001 * (pow(std::abs(MMM) / Mbar, 0.001)) * sign_MMM * (1.0 / Mbar);
				
				double denominator = df_dN * dN_dunp * pbProps_->mn_ + df_dV * dV_dusp * pbProps_->ms_ + df_dM * dM_dthetabp * pbProps_->mm_ + (df_dN * dN_dD + df_dV * dV_dD + df_dM * dM_dD) * (dD_dunp * pbProps_->mn_ + dD_dusp * pbProps_->ms_ + dD_dthetabp * pbProps_->mm_);
				delta_lamda_ = -1.0 * (pbProps_->Fyield_) / denominator;


				if ((N != 0.0) && (V == 0.0) && (M == 0.0)) {
					if (delta_lamda_<0.0) {

						delta_lamda_ = -1.0 * (pbProps_->Fyield_) / (df_dN * dN_dunp * pbProps_->mn_);
					}

				}


				
				if (delta_lamda_ < 0.0) {
					delta_lamda_ = prev_dellam_;
				}

				// we finish the calculation of the plastic multiplier, then we need to calculate the plastic deformation, damage variable and the new Fy(yield function value)
				int mark1 = 0; int mark2 = 0; int mark3 = 0;
				int signal1 = 0; int signal2 = 0; int signal3 = 0;
				double lamda0 = delta_lamda_;
				prev_npd_ = npd1_;
				prev_spd_ = spd1_;
				prev_bpr_ = bpr1_;
				npd1_ = normal_plastic_disp_;
				spd1_ = shear_plastic_disp_;
				bpr1_ = bending_plastic_ro_;

				// if the increment>0 
				if ((N != 0.0) && (pbProps_->mn_ != 0.0)) {
					if (N <= 0.0) {//tension let plastic deformation>0 also total disp>0
						//loading norm.x>0 or transx>0
						if (trans.x() > 0.0) {
							normal_plastic_disp_ = normal_plastic_disp_ + delta_lamda_ * std::abs(pbProps_->mn_);
							if (std::abs(delta_lamda_ * pbProps_->mn_) > dnorm) {
								normal_plastic_disp_ = npd1_ + dnorm;
								mark1 = 1;
							}
						}
						else {//unloading without new plastic deformation
							normal_plastic_disp_ = normal_plastic_disp_;
							signal1 = 1;
						}
					}
					else {//compression let plastic deformation<0 also total disp<0
						if (trans.x() <= 0.0) {//it means tensile loading
							normal_plastic_disp_ = normal_plastic_disp_ + delta_lamda_ * -1.0*(pbProps_->mn_);
							if (std::abs(delta_lamda_ * pbProps_->mn_) > dnorm) {
								normal_plastic_disp_ = npd1_ + trans.x();
								mark1 = 1;
							}
						}
						else {//unloading during tension without linitation
							normal_plastic_disp_ = normal_plastic_disp_;
								signal1 = 1;
						}
					}
				}
				if ((N != 0.0) && (pbProps_->mn_ == 0.0)) {
					normal_plastic_disp_ = normal_plastic_disp_ + (npd1_ - prev_npd_);
				}
				if (N == 0.0) {
					normal_plastic_disp_ = normal_plastic_disp_;
					signal1 = 1;
				}

				double increment_us = shear_total_disp_ - prev_std;
				// all shear total or plastic displacements>0 here, if increment_us<0 means unloading
			
				if ((V != 0.0) && (pbProps_->ms_ != 0.0)) {
					if (increment_us > 0.0) {//loading
						shear_plastic_disp_ = shear_plastic_disp_ + delta_lamda_ * pbProps_->ms_;
						if (std::abs(delta_lamda_ * pbProps_->ms_) > dtan) {
							shear_plastic_disp_ = spd1_ + dtan;
							mark2 = 1;
						}
					}
					else {//unloading
						shear_plastic_disp_ = shear_plastic_disp_ ;
						signal2 = 1;
					}
				}
				if ((V != 0.0) && (pbProps_->ms_ == 0.0)) {
					shear_plastic_disp_ = shear_plastic_disp_ + (spd1_ - prev_spd_);
				}
				if (V == 0.0) {
					shear_plastic_disp_ = shear_plastic_disp_;
					signal2 = 1;
				}
				
				double increment_btr = bending_total_ro_ - prev_btr;
				if ((M != 0.0) && (pbProps_->mm_ != 0.0)) {
					if (increment_btr>0.0) {
						bending_plastic_ro_ = bending_plastic_ro_ + delta_lamda_ * pbProps_->mm_;
						if (std::abs(delta_lamda_ * pbProps_->mm_) > dbang) {
							bending_plastic_ro_ = bpr1_ + dbang;
							mark3 = 1;
						}
					}
					else {
						bending_plastic_ro_ = bending_plastic_ro_ ;
						signal3 = 1;
					}
				}
				if ((M != 0.0) && (pbProps_->mm_ == 0.0)) {
					bending_plastic_ro_ = bending_plastic_ro_ + (bpr1_ - prev_bpr_);
				}
				if (M == 0.0) {
					bending_plastic_ro_ = bending_plastic_ro_;
					signal3 = 1;
				}
				

				
				prev_D_ = pbProps_->D_;
				prev_dellam_ = delta_lamda_;
				if (pbProps_->pb_F_.x() < 0.0) {//tension
					pbProps_->D_ = 1.0 - pow(2.718281828459, -(std::abs(normal_plastic_disp_) / pbProps_->ucn_ + std::abs(shear_plastic_disp_) / pbProps_->ucs_ + std::abs(bending_plastic_ro_) / pbProps_->ucb_));
				}
				
				else {
					pbProps_->D_ = 1.0 - pow(2.718281828459, -(std::abs(normal_plastic_disp_) / pbProps_->ucn1_ + std::abs(shear_plastic_disp_) / pbProps_->ucs_ + std::abs(bending_plastic_ro_) / pbProps_->ucb_));
				}
				

				//Save and Update contact strength properties
				pbProps_->pb_ten1_ = pbProps_->pb_ten_ * (1 - pbProps_->D_);
				pbProps_->pb_com1_ = pbProps_->pb_com_ * (1 - pbProps_->D_);
				pbProps_->pb_coh1_ = pbProps_->pb_coh_ * (1 - pbProps_->D_);
				//pbProps_->pb_mb1_ = pbProps_->pb_mb_ * (1 - pbProps_->D_);

				//Update kn, ks and Nbar Vbar Mbar
				pbProps_->pb_kn1_ = (1.0 - pbProps_->D_) * pbProps_->pb_kn_;
				pbProps_->pb_ks1_ = (1.0 - pbProps_->D_) * pbProps_->pb_ks_;
				if (N >= 0.0) {
					Nbar = pbProps_->pb_com1_ * pbArea;
				}
				else {
					Nbar = pbProps_->pb_ten1_ * pbArea;
				}
				Vbar = pbProps_->pb_coh1_ * pbArea;
				Mbar = pbProps_->pb_ten1_ * bi / br / br;

				//Update return forces and the bending moment, but we need to revise according to the mark1 2 3 value
				if (pbProps_->pb_F_.x() < 0.0) {
					if (signal1 == 0) {
						if (pbProps_->mn_ != 0.0) {
							if (mark1 == 0) {
								force_return_norm_ = (1.0 - pbProps_->D_) * delta_lamda_ * pbProps_->pb_kn_ * std::abs(pbProps_->mn_) * pbArea;
							}
							else {
								force_return_norm_ = (1.0 - pbProps_->D_) * dnorm * pbProps_->pb_kn_ * pbArea;
							}
						}
						else {
							force_return_norm_ = (1.0 - pbProps_->D_) * std::abs(npd1_ - prev_npd_) * pbProps_->pb_kn_ * pbArea;
						}
					}
					else {
						force_return_norm_ = 0.0;			
					}
				}
				else {
					if (signal1 == 0) {
						if (pbProps_->mn_ != 0.0) {
							if (mark1 == 0) {
								force_return_norm_ = -(1.0 - pbProps_->D_) * delta_lamda_ * pbProps_->pb_kn_ * std::abs(pbProps_->mn_) * pbArea;
							}
							else {
								force_return_norm_ = -(1.0 - pbProps_->D_) * dnorm * pbProps_->pb_kn_ * pbArea;
							}
						}
						else {
							force_return_norm_ = -(1.0 - pbProps_->D_) * std::abs(npd1_ - prev_npd_) * pbProps_->pb_kn_ * pbArea;
						}
					}
					else {
						force_return_norm_ = 0.0;
					}

				}

				

#ifdef THREED
				DVect return_forcex(1.0, 0.0, 0.0);//we have shear force in y and z direction, so we need to reassign shear force by ratio
				DVect return_forcey(0.0, 1.0, 0.0);
				DVect return_forcez(0.0, 0.0, 1.0);
				DVect return_momentx(1.0, 0.0, 0.0);
				DVect return_momenty(0.0, 1.0, 0.0);
				DVect return_momentz(0.0, 0.0, 1.0);
				if (pbProps_->pb_F_.y() != 0.0) {
					if (pbProps_->pb_F_.y() < 0.0) {
						if (pbProps_->ms_ != 0.0) {
							if (mark2 == 0) {
								force_return_tany_ = (1.0 - pbProps_->D_) * std::abs(bfs.ry()) / std::sqrt(bfs.ry() * bfs.ry() + bfs.rz() * bfs.rz()) * delta_lamda_ * pbProps_->pb_ks_ * pbProps_->ms_ * pbArea;
							}
							else {
								force_return_tany_ = (1.0 - pbProps_->D_) * std::abs(bfs.ry()) / std::sqrt(bfs.ry() * bfs.ry() + bfs.rz() * bfs.rz()) * std::abs(increment_us) * pbProps_->pb_ks_  * pbArea;
							}
						}
						else {
							force_return_tany_ = (1.0 - pbProps_->D_) * std::abs(bfs.ry()) / std::sqrt(bfs.ry() * bfs.ry() + bfs.rz() * bfs.rz()) * std::abs(spd1_ - prev_spd_) * pbProps_->pb_ks_  * pbArea;
						}
					}
					else { //Pfy>=0
						if (pbProps_->ms_ != 0.0) {
							if (mark2 == 0) {
								force_return_tany_ = -(1.0 - pbProps_->D_) * std::abs(bfs.ry()) / std::sqrt(bfs.ry() * bfs.ry() + bfs.rz() * bfs.rz()) * delta_lamda_ * pbProps_->pb_ks_ * pbProps_->ms_ * pbArea;
							}
							else {
								force_return_tany_ = -(1.0 - pbProps_->D_) * std::abs(bfs.ry()) / std::sqrt(bfs.ry() * bfs.ry() + bfs.rz() * bfs.rz()) * std::abs(increment_us) * pbProps_->pb_ks_ * pbArea;

							}
						}
						else {
							force_return_tany_ = -(1.0 - pbProps_->D_) * std::abs(bfs.ry()) / std::sqrt(bfs.ry() * bfs.ry() + bfs.rz() * bfs.rz()) * std::abs(spd1_ - prev_spd_) * pbProps_->pb_ks_  * pbArea;
						}
					}
				}

				if (pbProps_->pb_F_.z() != 0.0) {
					if (pbProps_->pb_F_.z() < 0.0) {
						if (pbProps_->ms_ != 0.0) {
							if (mark2 == 0) {
								force_return_tanz_ = (1.0 - pbProps_->D_) * std::abs(bfs.rz()) / std::sqrt(bfs.ry() * bfs.ry() + bfs.rz() * bfs.rz()) * delta_lamda_ * pbProps_->pb_ks_ * pbProps_->ms_ * pbArea;
							}
							else {
								force_return_tanz_ = (1.0 - pbProps_->D_) * std::abs(bfs.rz()) / std::sqrt(bfs.ry() * bfs.ry() + bfs.rz() * bfs.rz()) * std::abs(increment_us) * pbProps_->pb_ks_ * pbArea;
							}
						}
						else {
							force_return_tanz_ = (1.0 - pbProps_->D_) * std::abs(bfs.rz()) / std::sqrt(bfs.ry() * bfs.ry() + bfs.rz() * bfs.rz()) * std::abs(spd1_ - prev_spd_) * pbProps_->pb_ks_ * pbArea;
						}
					}
					else {
						if (pbProps_->ms_ != 0.0) {
							if (mark2 == 0) {
								force_return_tanz_ = -(1.0 - pbProps_->D_) * std::abs(bfs.rz()) / std::sqrt(bfs.ry() * bfs.ry() + bfs.rz() * bfs.rz()) * delta_lamda_ * pbProps_->pb_ks_ * pbProps_->ms_ * pbArea;
							}
							else {
								force_return_tanz_ = -(1.0 - pbProps_->D_) * std::abs(bfs.rz()) / std::sqrt(bfs.ry() * bfs.ry() + bfs.rz() * bfs.rz()) * std::abs(increment_us) * pbProps_->pb_ks_ * pbArea;
							}
						}
						else {
							force_return_tanz_ = -(1.0 - pbProps_->D_) * std::abs(bfs.rz()) / std::sqrt(bfs.ry() * bfs.ry() + bfs.rz() * bfs.rz()) * std::abs(spd1_ - prev_spd_) * pbProps_->pb_ks_ * pbArea;
						}
					}
				}
				if (signal2 == 1) {
					force_return_tanz_ = 0.0;
					force_return_tany_ = 0.0;
				}
				if (pbProps_->pb_M_.y() != 0.0) {
					if (pbProps_->pb_M_.y() < 0.0) {
						if (pbProps_->mm_ != 0.0) {
							if (mark3 == 0) {
								moment_return_by_ = (1.0 - pbProps_->D_) * std::abs(bm.ry()) / std::sqrt(bm.ry() * bm.ry() + bm.rz() * bm.rz()) * delta_lamda_ * pbProps_->pb_kn_ * pbProps_->mm_ * bi;
							}
							else {
								moment_return_by_ = (1.0 - pbProps_->D_) * std::abs(bm.ry()) / std::sqrt(bm.ry() * bm.ry() + bm.rz() * bm.rz()) * std::abs(increment_btr) * pbProps_->pb_kn_ * bi;
							}
						}
						else {
							moment_return_by_ = (1.0 - pbProps_->D_) * std::abs(bm.ry()) / std::sqrt(bm.ry() * bm.ry() + bm.rz() * bm.rz()) * std::abs(bpr1_ - prev_bpr_) * pbProps_->pb_kn_ * bi;
						}	
					}
					else {
						if (pbProps_->mm_ != 0.0) {
							if (mark3 == 0) {
								moment_return_by_ = -(1.0 - pbProps_->D_) * std::abs(bm.ry()) / std::sqrt(bm.ry() * bm.ry() + bm.rz() * bm.rz()) * delta_lamda_ * pbProps_->pb_kn_ * pbProps_->mm_ * bi;
							}
							else {
								moment_return_by_ = -(1.0 - pbProps_->D_) * std::abs(bm.ry()) / std::sqrt(bm.ry() * bm.ry() + bm.rz() * bm.rz()) * std::abs(increment_btr) * pbProps_->pb_kn_ * bi;
							}
						}
						else {
							moment_return_by_ = -(1.0 - pbProps_->D_) * std::abs(bm.ry()) / std::sqrt(bm.ry() * bm.ry() + bm.rz() * bm.rz()) * std::abs(bpr1_ - prev_bpr_) * pbProps_->pb_kn_ * bi;
						}
					}
				}
				if (pbProps_->pb_M_.z() != 0.0) {
					if (pbProps_->pb_M_.z() < 0.0) {
						if (pbProps_->mm_ != 0.0) {
							if (mark3 == 0) {
								moment_return_bz_ = (1.0 - pbProps_->D_) * std::abs(bm.rz()) / std::sqrt(bm.ry() * bm.ry() + bm.rz() * bm.rz()) * delta_lamda_ * pbProps_->pb_kn_ * pbProps_->mm_ * bi;
							}
							else {
								moment_return_bz_ = (1.0 - pbProps_->D_) * std::abs(bm.rz()) / std::sqrt(bm.ry() * bm.ry() + bm.rz() * bm.rz()) * std::abs(increment_btr) * pbProps_->pb_kn_ * bi;
							}
						}
						else {
							moment_return_bz_ = (1.0 - pbProps_->D_) * std::abs(bm.rz()) / std::sqrt(bm.ry() * bm.ry() + bm.rz() * bm.rz()) * std::abs(bpr1_ - prev_bpr_) * pbProps_->pb_kn_ * bi;
						}			
					}
					else {
						if (pbProps_->mm_ != 0.0) {
							if (mark3 == 0) {
								moment_return_bz_ = -(1.0 - pbProps_->D_) * std::abs(bm.rz()) / std::sqrt(bm.ry() * bm.ry() + bm.rz() * bm.rz()) * delta_lamda_ * pbProps_->pb_kn_ * pbProps_->mm_ * bi;
							}
							else {
								moment_return_bz_ = -(1.0 - pbProps_->D_) * std::abs(bm.rz()) / std::sqrt(bm.ry() * bm.ry() + bm.rz() * bm.rz()) * std::abs(increment_btr) * pbProps_->pb_kn_ * bi;
							}
						}
						else {
							moment_return_bz_ = -(1.0 - pbProps_->D_) * std::abs(bm.rz()) / std::sqrt(bm.ry() * bm.ry() + bm.rz() * bm.rz()) * std::abs(bpr1_ - prev_bpr_) * pbProps_->pb_kn_ * bi;
						}
					}
				}

				if (signal3 == 1) {
					moment_return_bz_ = 0.0;
					moment_return_by_ = 0.0;
				}
				//Return force and moment vector 
				//update or revise the current bond force and moment, you need to use the return revised force calculated above , here 3 directions are all revised
				/*DAVect Mx(pbProps_->pb_M_.rx(), 0.0, 0.0);
				double twist_moment;
				if (pbProps_->pb_M_.rx() < 0.0) {
					twist_moment = Mx.mag();
				}
				else {
					twist_moment = -1.0*Mx.mag();
				}*/
				
				pbProps_->pb_F_ = pbProps_->pb_F_ * (1.0 - pbProps_->D_) / (1.0 - prev_D_) + return_forcex * force_return_norm_ + return_forcey * force_return_tany_ + return_forcez * force_return_tanz_;
				pbProps_->pb_M_ = pbProps_->pb_M_ * (1.0 - pbProps_->D_) / (1.0 - prev_D_) + return_momenty * moment_return_by_ + return_momentz * moment_return_bz_;
				double N1 = pbProps_->pb_F_.rx();
				DVect temp1(pbProps_->pb_F_);
				DAVect temp2(pbProps_->pb_M_);
				temp1.rx() = 0.0;
				temp2.rx() = 0.0;
				double V1 = temp1.mag();
				double M1 = temp2.mag();
				BendingMo_ = M1;
				BendingMotilde = BendingMo_ / br;
				double Nbar1;
				if (pbProps_->pb_F_.x() < 0.0) {
					Nbar1 = pbProps_->pb_ten_ * (1.0 - prev_D_) * pbArea;
				}
				else {
					Nbar1 = pbProps_->pb_com_ * (1.0 - prev_D_) * pbArea;
				}
				double Vbar1 = pbProps_->pb_coh_ * (1.0 - prev_D_) * pbArea;
				double Mbar1 = pbProps_->pb_ten_ * (1.0 - prev_D_) * bi / br / br;
				Vy_ = pbProps_->pb_F_.y();
				Vz_ = pbProps_->pb_F_.z();
				BendingMoy_ = pbProps_->pb_M_.y();
				BendingMoz_ = pbProps_->pb_M_.z();
				pbProps_->Fy_ = pow(std::abs(BendingMotilde / Mbar1), 1.001) + pow((N1 / Nbar1), 2.0) + (pow((V1 / Vbar1), 4.0) / (1.0 - pow((N1 / Nbar1), 2.0))) - 1.0;
			
				//calculate the yield value again to judge if it needs to revise again,
				//double Ntemp = pbProps_->pb_F_.x();
				//double Vtemp = 

#else
				DVect return_forcex(1.0, 0.0);//Return force vector 2D case, !!!!!!!!I havent revise this part!!!!!!!!!!!!!!
				DVect return_forcey(0.0, 1.0);
				if (pbProps_->pb_F_.y() < 0.0) {
					force_return_tany_ = (1.0 - pbProps_->D_) * delta_lamda_ * pbProps_->pb_ks_ * pbProps_->ms_ * pbArea;
				}
				else {
					force_return_tany_ = -(1.0 - pbProps_->D_) * delta_lamda_ * pbProps_->pb_ks_ * pbProps_->ms_ * pbArea;
				}
				pbProps_->pb_F_ = pbhProps_->pb_F_ * (1.0 - pbProps_->D_) / (1.0 - prev_D_) + return_forcex * force_return_norm_ + return_forcey * force_return_tany_;
#endif
			} //Fyield > 0


			

			/* check bond failure */
			if (state->canFail_) {
				/* maximum stresses */
				double dbend = sqrt(pbProps_->pb_M_.y()*pbProps_->pb_M_.y() + pbProps_->pb_M_.z()*pbProps_->pb_M_.z());
				double dtwist = pbProps_->pb_M_.x();
				DVect bfs(pbProps_->pb_F_);
				bfs.rx() = 0.0;
				double dbfs = bfs.mag();
				double nsmax = -(pbProps_->pb_F_.x() / pbArea);  //+  pbProps_->pb_mcf_ * dbend * br/bi; 
				double ssmax = dbfs / pbArea;  //+  pbProps_->pb_mcf_ * std::abs(dtwist) * 0.5* br/bi;
				double ss;


				if (pbProps_->D_ > 0.99 ) {
					double se = pbStrainEnergy(); // bond strain energy at the onset of failure
					pbProps_->pb_state_ = 1;
                    pbProps_->pb_F_.fill(0.0);
                    pbProps_->pb_M_.fill(0.0);
                    if (cmEvents_[fBondBreak] >= 0) {
						auto c = state->getContact();
						//QVariant p1;
						//IContact* cp = const_cast<IContact*>(state->getContact());
						//TPtr<IThing> t(cp->getIThing());
						//p1.setValue(t);
						std::vector<fish::Parameter> arg = { fish::Parameter(c->getIThing()),
							                                 //fish::Parameter((qint64)abc),
															 fish::Parameter((qint64)pbProps_->pb_state_),
															 fish::Parameter(pbProps_->pb_ten1_),
															 fish::Parameter(se) };
						//auto c = state->getContact();
                        IFishCallList *fi = const_cast<IFishCallList*>(state->getProgram()->findInterface<IFishCallList>());
                        fi->setCMFishCallArguments(c,arg,cmEvents_[fBondBreak]);   //call fish function to break the bond
                    }
                }
            }
        }// pb_F calculation
		

        // Compute energies
        if (state->trackEnergy_) {
            assert(energies_);
            energies_->estrain_ = 0.0;
            energies_->epbstrain_ = 0.0;
            if (kn_)
                energies_->estrain_ = 0.5 * lin_F_.x()* lin_F_.x() /kn_;
            if (ks_) {
                DVect s = lin_F_;
                s.rx() = 0.0;
                double smag2 = s.mag2();
                energies_->estrain_ += 0.5* smag2 / ks_ ;

                if (lin_S_) {
                    ss_F_old.rx() = 0.0;
                    DVect avg_F_s = (s + ss_F_old)*0.5;
                    DVect u_s_el =  (s - ss_F_old) / ks_;
                    energies_->eslip_ -= std::min(0.0,(avg_F_s | (u_s + u_s_el)));
                }
            }
            if (dpProps_)
                energies_->edashpot_ -= dpProps_->dp_F_ | trans;
            if (doPb)
                energies_->epbstrain_ = pbStrainEnergy();
        }
		//}


		//state->getMechanicalContact()->updateResultingTorquesLocal(state->force_, &state->momentOn1_, &state->momentOn2_);
        assert(lin_F_ == lin_F_);
        return checkActivity(state->gap_);
    }
    bool ContactModelDBondrotation::thermalCoupling(ContactModelMechanicalState *ms,ContactModelThermalState *ts, IContactThermal *ct,const double &) {
        bool ret = false;
        if (!pbProps_) return ret;
        if (pbProps_->pb_state_ < 4) return ret;
        int idx = ct->getModel()->getContactModel()->isProperty("thexp");
        if (idx<=0 ) return ret;

        double thexp = (ct->getModel()->getContactModel()->getProperty(idx)).toDouble();
        double length = ts->length_;
        double delTemp =ts->tempInc_;
        double delUn = length * thexp * delTemp;
        if (delUn == 0.0) return ret;

        double dthick = 0.0;
        double Cmin1 = ms->end1Curvature_.x();
        double Cmax1 = ms->end1Curvature_.y();
        double Cmin2 = ms->end2Curvature_.x();
        double Cmax2 = ms->end2Curvature_.y();

        Cmin2;
        if (Cmin1 == 0.0)
            //dthick = ms->end1Extent_.x();
		    dthick = 1.0;

        double br = pbProps_->pb_rmul_ * 1.0 / std::max(Cmax1,Cmax2);
        if (userArea_)
#ifdef THREED
            br = std::sqrt(userArea_ / dPi);
#else
            br = userArea_ / 2.0;
#endif
        double br2 = br*br;
        double pbArea = dthick <= 0.0 ? dPi*br2 : 2.0*br*dthick;
        //
        DVect finc(0.0);
        finc.rx() = 1.0;
        finc *= pbProps_->pb_kn1_*pbArea*delUn;
        pbProps_->pb_F_ += finc;

        //ms->force_ += finc;

        // The state force has been updated - update the state with the resulting torques
        //ms->getMechanicalContact()->updateResultingTorquesLocal(ms->force_,&ms->momentOn1_,&ms->momentOn2_);

        ret = true;
        return ret;
    }

    void ContactModelDBondrotation::setForce(const DVect &v,IContact *c) {
        lin_F(v);
        if (v.x() > 0)
            rgap_ = c->getGap() + v.x() / kn_;
    }

    void ContactModelDBondrotation::propagateStateInformation(IContactModelMechanical* old,const CAxes &oldSystem,const CAxes &newSystem) {
        // Only do something if the contact model is of the same type
        if (old->getContactModel()->getName().compare("dbondrotation",Qt::CaseInsensitive) == 0 && !isBonded()) {
            ContactModelDBondrotation *oldCm = (ContactModelDBondrotation *)old;
#ifdef THREED
            // Need to rotate just the shear component from oldSystem to newSystem

            // Step 1 - rotate oldSystem so that the normal is the same as the normal of newSystem
            DVect axis = oldSystem.e1() & newSystem.e1();
            double c, ang, s;
            DVect re2;
            if (!checktol(axis.abs().maxComp(),0.0,1.0,1000)) {
                axis = axis.unit();
                c = oldSystem.e1()|newSystem.e1();
                if (c > 0)
                    c = std::min(c,1.0);
                else
                    c = std::max(c,-1.0);
                ang = acos(c);
                s = sin(ang);
                double t = 1. - c;
                DMatrix<3,3> rm;
                rm.get(0,0) = t*axis.x()*axis.x() + c;
                rm.get(0,1) = t*axis.x()*axis.y() - axis.z()*s;
                rm.get(0,2) = t*axis.x()*axis.z() + axis.y()*s;
                rm.get(1,0) = t*axis.x()*axis.y() + axis.z()*s;
                rm.get(1,1) = t*axis.y()*axis.y() + c;
                rm.get(1,2) = t*axis.y()*axis.z() - axis.x()*s;
                rm.get(2,0) = t*axis.x()*axis.z() - axis.y()*s;
                rm.get(2,1) = t*axis.y()*axis.z() + axis.x()*s;
                rm.get(2,2) = t*axis.z()*axis.z() + c;
                re2 = rm*oldSystem.e2();
            }
            else
                re2 = oldSystem.e2();
            // Step 2 - get the angle between the oldSystem rotated shear and newSystem shear
            axis = re2 & newSystem.e2();
            DVect2 tpf;
            DMatrix<2,2> m;
            if (!checktol(axis.abs().maxComp(),0.0,1.0,1000)) {
                axis = axis.unit();
                c = re2|newSystem.e2();
                if (c > 0)
                    c = std::min(c,1.0);
                else
                    c = std::max(c,-1.0);
                ang = acos(c);
                if (!checktol(axis.x(),newSystem.e1().x(),1.0,100))
                    ang *= -1;
                s = sin(ang);
                m.get(0,0) = c;
                m.get(1,0) = s;
                m.get(0,1) = -m.get(1,0);
                m.get(1,1) = m.get(0,0);
                tpf = m*DVect2(oldCm->lin_F_.y(),oldCm->lin_F_.z());
            } else {
                m.get(0,0) = 1.;
                m.get(0,1) = 0.;
                m.get(1,0) = 0.;
                m.get(1,1) = 1.;
                tpf = DVect2(oldCm->lin_F_.y(),oldCm->lin_F_.z());
            }
            DVect pforce = DVect(0,tpf.x(),tpf.y());
#else
            oldSystem;
            newSystem;
            DVect pforce = DVect(0,oldCm->lin_F_.y());
#endif
            for (int i=1; i<dim; ++i)
                lin_F_.rdof(i) += pforce.dof(i);
			if (lin_mode_ && oldCm->lin_mode_)
				lin_F_.rx() = oldCm->lin_F_.x();
            oldCm->lin_F_ = DVect(0.0);
            if (dpProps_ && oldCm->dpProps_) {
#ifdef THREED
                tpf = m*DVect2(oldCm->dpProps_->dp_F_.y(),oldCm->dpProps_->dp_F_.z());
                pforce = DVect(oldCm->dpProps_->dp_F_.x(),tpf.x(),tpf.y());
#else
                pforce = oldCm->dpProps_->dp_F_;
#endif
                dpProps_->dp_F_ += pforce;
                oldCm->dpProps_->dp_F_ = DVect(0.0);
            }
            if(oldCm->getEnergyActivated()) {
                activateEnergy();
                energies_->estrain_ = oldCm->energies_->estrain_;
                energies_->eslip_ = oldCm->energies_->eslip_;
                energies_->edashpot_ = oldCm->energies_->edashpot_;
                energies_->epbstrain_ = oldCm->energies_->epbstrain_;
                oldCm->energies_->estrain_ = 0.0;
                oldCm->energies_->edashpot_ = 0.0;
                oldCm->energies_->eslip_ = 0.0;
                oldCm->energies_->epbstrain_ = 0.0;
            }
            rgap_ = oldCm->rgap_;
        }
        assert(lin_F_ == lin_F_);
    }

    void ContactModelDBondrotation::setNonForcePropsFrom(IContactModel *old) {
        // Only do something if the contact model is of the same type
        if (old->getName().compare("dbondrotation",Qt::CaseInsensitive) == 0 && !isBonded()) {
            ContactModelDBondrotation *oldCm = (ContactModelDBondrotation *)old;
            kn_ = oldCm->kn_;
            ks_ = oldCm->ks_;
            fric_ = oldCm->fric_;
            lin_mode_ = oldCm->lin_mode_;
            rgap_ = oldCm->rgap_;
            userArea_ = oldCm->userArea_;

            if (oldCm->dpProps_) {
                if (!dpProps_)
                    dpProps_ = NEWC(dpProps());
                dpProps_->dp_nratio_ = oldCm->dpProps_->dp_nratio_;
                dpProps_->dp_sratio_ = oldCm->dpProps_->dp_sratio_;
                dpProps_->dp_mode_ = oldCm->dpProps_->dp_mode_;
            }

            if (oldCm->pbProps_) {
                if (!pbProps_)
                    pbProps_ = NEWC(pbProps());
                pbProps_->pb_rmul_ = oldCm->pbProps_->pb_rmul_;
                pbProps_->pb_kn1_ = oldCm->pbProps_->pb_kn1_;
                pbProps_->pb_ks1_ = oldCm->pbProps_->pb_ks1_;
                pbProps_->pb_mcf_ = oldCm->pbProps_->pb_mcf_;
                //pbProps_->pb_fa_ = oldCm->pbProps_->pb_fa_;
                pbProps_->pb_state_ = oldCm->pbProps_->pb_state_;
                pbProps_->pb_coh1_ = oldCm->pbProps_->pb_coh1_;
                pbProps_->pb_ten1_ = oldCm->pbProps_->pb_ten1_;
				pbProps_->pb_com1_ = oldCm->pbProps_->pb_com1_;
				//pbProps_->pb_mb1_ = oldCm->pbProps_->pb_mb1_;
				pbProps_->ucn_ = oldCm->pbProps_->ucn_;
				pbProps_->ucn1_ = oldCm->pbProps_->ucn1_;
				pbProps_->ucs_ = oldCm->pbProps_->ucs_;
				pbProps_->ucb_ = oldCm->pbProps_->ucb_;
				//pbProps_->alpha_ = oldCm->pbProps_->alpha_;
				//pbProps_->timing_ = oldCm->pbProps_->timing_;
				pbProps_->D_ = oldCm->pbProps_->D_;
				pbProps_->Fyield_ = oldCm->pbProps_->Fyield_;
				pbProps_->Fy_ = oldCm->pbProps_->Fy_;
                pbProps_->pbTransStiff_ = oldCm->pbProps_->pbTransStiff_;
                pbProps_->pbAngStiff_ = oldCm->pbProps_->pbAngStiff_;
            }
        }
    }

    DVect ContactModelDBondrotation::getForce(const IContactMechanical *) const {
        DVect ret(lin_F_);
        if (dpProps_)
            ret += dpProps_->dp_F_;
		if (pbProps_)

            ret += pbProps_->pb_F_;
        return ret;
    }

    DAVect ContactModelDBondrotation::getMomentOn1(const IContactMechanical *c) const {
        DVect force = getForce(c);
        DAVect ret(0.0);
        if (pbProps_)
            ret = pbProps_->pb_M_;
        c->updateResultingTorqueOn1Local(force,&ret);
        return ret;
    }

    DAVect ContactModelDBondrotation::getMomentOn2(const IContactMechanical *c) const {
        DVect force = getForce(c);
        DAVect ret(0.0);
        if (pbProps_)
            ret = pbProps_->pb_M_;
        c->updateResultingTorqueOn2Local(force,&ret);
        return ret;
    }

    DVect3 ContactModelDBondrotation::pbData(const IContactMechanical *c) const {
        //double Cmin1 = c->getEnd1Curvature().x();
        double Cmax1 = c->getEnd1Curvature().y();
        double Cmax2 = c->getEnd2Curvature().y();
        //double dthick = (Cmin1 == 0.0) ? c->getEnd1Extent().x() : 0.0;
        double br = pbProps_->pb_rmul_ * 1.0 / std::max(Cmax1,Cmax2);
        if (userArea_)
#ifdef THREED
            br = std::sqrt(userArea_ / dPi);
#else
            br = userArea_ / 2.0;
#endif
        double br2 = br*br;
#ifdef TWOD
		double pbArea = 2.0*br;
		double bi = 2.0*br*br2 / 3.0;
#else
		double pbArea = dPi * br2;
		double bi = 0.25*pbArea*br2;
#endif
        return DVect3(pbArea,bi,br);
    }

    DVect2 ContactModelDBondrotation::pbSMax(const IContactMechanical *c) const {
        DVect3 data = pbData(c);
        double pbArea = data.x();
        double bi = data.y();
        double br = data.z();
        /* maximum stresses */
        double dbend = sqrt(pbProps_->pb_M_.y()*pbProps_->pb_M_.y() + pbProps_->pb_M_.z()*pbProps_->pb_M_.z());
        double dtwist = pbProps_->pb_M_.x();
        DVect bfs(pbProps_->pb_F_);
        bfs.rx() = 0.0;
        double dbfs = bfs.mag();
		double nsmax = -(pbProps_->pb_F_.x() / pbArea);  //+  pbProps_->pb_mcf_ * dbend * br/bi;
		double ssmax = dbfs / pbArea; // +  pbProps_->pb_mcf_ * std::abs(dtwist) * 0.5* br/bi;
		//double nor_ef = -(pbProps_->pb_F_.x() / pbArea);  //+  pbProps_->pb_mcf_ * dbend * br/bi;
		//double she_ef = dbfs / pbArea;
        return DVect2(nsmax,ssmax);
		//return DVect2(nor_ef,she_ef);
    }

    /*double ContactModelDBondrotation::pbShearStrength(const double &pbArea) const {
        if (!pbProps_) return 0.0;
        double sig = -1.0*pbProps_->pb_F_.x() / pbArea;
        double nstr1 = pbProps_->pb_state_ > 3 ? pbProps_->pb_ten1_ : 0.0;
		double nstr2 = pbProps_->pb_state_ > 3 ? (-1.0)*pbProps_->pb_com1_ : 0.0;
		if (sig >= 0.0) {
			return sig <= nstr1 ? pbProps_->pb_coh1_ - std::tan(dDegrad*pbProps_->pb_fa_)*sig
				: pbProps_->pb_coh1_ - std::tan(dDegrad*pbProps_->pb_fa_)*nstr1;
		}
		if (sig < 0.0) {
			return sig >= nstr2 ? pbProps_->pb_coh1_ - std::tan(dDegrad*pbProps_->pb_fa_)*sig
				: pbProps_->pb_coh1_ - std::tan(dDegrad*pbProps_->pb_fa_)*nstr2;
		}
    }*/

    void ContactModelDBondrotation::setDampCoefficients(const double &mass,double *vcn,double *vcs) {
        *vcn = dpProps_->dp_nratio_ * 2.0 * sqrt(mass*(kn_));
        *vcs = dpProps_->dp_sratio_ * 2.0 * sqrt(mass*(ks_));
    }

} // namespace itascaxd
// EoF
