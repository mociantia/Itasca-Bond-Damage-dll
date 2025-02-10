#pragma once
// contactmodeldbondrotation.h

#include "contactmodel/src/contactmodelmechanical.h"

#ifdef DBondrotation_LIB
#  define DBondrotation_EXPORT EXPORT_TAG
#elif defined(NO_MODEL_IMPORT)
#  define DBondrotation_EXPORT
#else
#  define DBondrotation_EXPORT IMPORT_TAG
#endif

namespace cmodelsxd {
    using namespace itasca;

    class ContactModelDBondrotation : public ContactModelMechanical {
    public:
        DBondrotation_EXPORT ContactModelDBondrotation();
        DBondrotation_EXPORT virtual ~ContactModelDBondrotation();
        virtual void                     copy(const ContactModel *c) override;
        virtual void                     archive(ArchiveStream &);
        virtual QString  getName() const { return "dbondrotation"; }
        virtual void     setIndex(int i) { index_=i;}
        virtual int      getIndex() const {return index_;}

		enum PropertyKeys {
			kwLinKn = 1
			, kwLinKs
			, kwLinFric
			, kwLinF
			, kwLinS
			, kwLinMode
			, kwRGap
			, kwEmod
			, kwKRatio
			, kwDpNRatio
			, kwDpSRatio
			, kwDpMode
			, kwDpF
			, kwPbState
			, kwPbRMul
			, kwPbKn
			, kwPbKs
			, kwPbKn1
			, kwPbKs1
			, kwPbMcf
			, kwPbTStrength
			, kwPbCStrength
            //, kwPbMBStrength
			, kwPbTStrength1
			, kwPbCStrength1
			//, kwPbMBStrength1
            , kwPbCoh
			, kwPbCoh1
            //, kwPbFa
            //, kwPbSig
            //, kwPbTau
            , kwPbF
            , kwPbM
            , kwPbRadius
            , kwPbEmod
            , kwPbKRatio
            , kwUserArea
            , kwMN
            , kwMS
			, kwMM
			, kwFYield
			, kwUCN
			, kwUCN1
			, kwUCS
			, kwUCB
			, kwDamagevar
			, kwFy
			, kwNtd
			, kwStdy
			, kwStdz
			, kwStd
			, kwBtry
			, kwBtrz
			, kwBtr
			, kwNpd
			, kwSpd
			, kwBpr
			, kwBendingMy
			, kwBendingMz
			, kwBendingM
			, kwSheary
			, kwShearz
			//, kwALPHA
			//, kwTIMING
        };

        virtual QString  getProperties() const {
			return "kn"
				",ks"
				",fric"
				",lin_force"
				",lin_slip"
				",lin_mode"
				",rgap"
				",emod"
				",kratio"
				",dp_nratio"
				",dp_sratio"
				",dp_mode"
				",dp_force"
				",pb_state"
				",pb_rmul"
				",pb_kn"
				",pb_ks"
				",pb_kn1"
				",pb_ks1"
				",pb_mcf"
				",pb_ten"
				",pb_com"
				//",pb_mb"
				",pb_ten1"
				",pb_com1"
				//",pb_mb1"
				",pb_coh"
				",pb_coh1"
				//",pb_fa"
				//",pb_sigma"
				//",pb_tau"
				",pb_force"
				",pb_moment"
				",pb_radius"
				",pb_emod"
				",pb_kratio"
				",user_area"
				",mn"
				",ms"
				",mm"
				",Fyield"
				",ucn"
				",ucn1"
				",ucs"
				",ucb"
				",D_var"
				",Fy"
				",Ntd"
				",Stdy"
				",Stdz"
				",Std"
				",Btry"
				",Btrz"
				",Btr"
				",Npd"
				",Spd"
				",Bpr"
				",BendingMy"
				",BendingMz"
				",BendingM"
				",Sheary"
				",Shearz";
				//",alpha"
				//",timing";
        }

        enum EnergyKeys { kwEStrain=1,kwESlip,kwEDashpot,kwEPbStrain};
        virtual QString  getEnergies() const { return "energy-strain,energy-slip,energy-dashpot,energy-pbstrain";}
        virtual double   getEnergy(uint i) const;  // Base 1
        virtual bool     getEnergyAccumulate(uint i) const; // Base 1
        virtual void     setEnergy(uint i,const double &d); // Base 1
        virtual void     activateEnergy() { if (energies_) return; energies_ = NEWC(Energies());}
        virtual bool     getEnergyActivated() const {return (energies_ !=0);}

        enum FishCallEvents {fActivated=0,fBondBreak,fSlipChange};
        virtual QString  getFishCallEvents() const { return "contact_activated,bond_break,slip_change"; }
        virtual QVariant getProperty(uint i,const IContact *con=0) const;
        virtual bool     getPropertyGlobal(uint i) const;
        virtual bool     setProperty(uint i,const QVariant &v,IContact *con=0);
        virtual bool     getPropertyReadOnly(uint i) const;

        virtual bool     supportsInheritance(uint i) const;
        virtual bool     getInheritance(uint i) const { assert(i<32); quint32 mask = to<quint32>(1 << i);  return (inheritanceField_ & mask) ? true : false; }
        virtual void     setInheritance(uint i,bool b) { assert(i<32); quint32 mask = to<quint32>(1 << i);  if (b) inheritanceField_ |= mask;  else inheritanceField_ &= ~mask; }

        enum MethodKeys { kwDeformability=1
                        , kwPbDeformability
                        , kwPbBond
                        , kwPbUnbond
                        , kwArea
        };

        virtual QString  getMethods() const {
            return "deformability"
                   ",pb_deformability"
                   ",bond"
                   ",unbond"
                   ",area";
        }

        virtual QString  getMethodArguments(uint i) const;

        virtual bool     setMethod(uint i,const QVector<QVariant> &vl,IContact *con=0); // Base 1 - returns true if timestep contributions need to be updated

        virtual uint     getMinorVersion() const;

        virtual bool    validate(ContactModelMechanicalState *state,const double &timestep);
        virtual bool    endPropertyUpdated(const QString &name,const IContactMechanical *c);
        virtual bool    forceDisplacementLaw(ContactModelMechanicalState *state,const double &timestep);
        virtual DVect2  getEffectiveTranslationalStiffness() const { DVect2 ret = effectiveTranslationalStiffness_; if(pbProps_) ret+= pbProps_->pbTransStiff_ ;return ret;}
        virtual DAVect  getEffectiveRotationalStiffness() const {if (!pbProps_) return DAVect(0.0); return pbProps_->pbAngStiff_;}

        virtual bool thermalCoupling(ContactModelMechanicalState *, ContactModelThermalState * , IContactThermal *,const double &);

        virtual ContactModelDBondrotation *clone() const override { return NEWC(ContactModelDBondrotation()); }
        virtual double              getActivityDistance() const {return rgap_;}
        virtual bool                isOKToDelete() const { return !isBonded(); }
        virtual void                resetForcesAndMoments() { lin_F(DVect(0.0)); dp_F(DVect(0.0)); pbF(DVect(0.0)); pbM(DAVect(0.0));  if (energies_) { energies_->estrain_ = 0.0;  if (energies_) energies_->epbstrain_ = 0.0;}}
        virtual void                setForce(const DVect &v,IContact *c);
        virtual void                setArea(const double &d) { userArea_ = d; }
		virtual double              getArea() const { return userArea_; }

        virtual bool     checkActivity(const double &gap) { return (gap <= rgap_ || isBonded()); }

        virtual bool     isSliding() const { return lin_S_; }
        virtual bool     isBonded() const { return pbProps_ ? (pbProps_->pb_state_==4 ) : false; }
		virtual void     unbond() { if (pbProps_) pbProps_->pb_state_ = 0; }
        virtual void     propagateStateInformation(IContactModelMechanical* oldCm,const CAxes &oldSystem=CAxes(),const CAxes &newSystem=CAxes());
        virtual void     setNonForcePropsFrom(IContactModel *oldCM);

        const double & kn() const {return kn_;}
        void           kn(const double &d) {kn_=d;}
        const double & ks() const {return ks_;}
        void           ks(const double &d) {ks_=d;}
        const double & fric() const {return fric_;}
        void           fric(const double &d) {fric_=d;}
        const DVect &  lin_F() const {return lin_F_;}
        void           lin_F(const DVect &f) { lin_F_=f;}
        bool           lin_S() const {return lin_S_;}
        void           lin_S(bool b) { lin_S_=b;}
        uint           lin_mode() const {return lin_mode_;}
        void           lin_mode(uint i) { lin_mode_=i;}
        const double & rgap() const {return rgap_;}
        void           rgap(const double &d) {rgap_=d;}

        bool     hasDamping() const {return dpProps_ ? true : false;}
        double   dp_nratio() const {return (hasDamping() ? (dpProps_->dp_nratio_) : 0.0);}
        void     dp_nratio(const double &d) { if(!hasDamping()) return; dpProps_->dp_nratio_=d;}
        double   dp_sratio() const {return hasDamping() ? dpProps_->dp_sratio_: 0.0;}
        void     dp_sratio(const double &d) { if(!hasDamping()) return; dpProps_->dp_sratio_=d;}
        int      dp_mode() const {return hasDamping() ? dpProps_->dp_mode_: -1;}
        void     dp_mode(int i) { if(!hasDamping()) return; dpProps_->dp_mode_=i;}
        DVect    dp_F() const {return hasDamping() ? dpProps_->dp_F_: DVect(0.0);}
        void     dp_F(const DVect &f) { if(!hasDamping()) return; dpProps_->dp_F_=f;}

        bool    hasEnergies() const {return energies_ ? true:false;}
        double  estrain() const {return hasEnergies() ? energies_->estrain_: 0.0;}
        void    estrain(const double &d) { if(!hasEnergies()) return; energies_->estrain_=d;}
        double  eslip() const {return hasEnergies() ? energies_->eslip_: 0.0;}
        void    eslip(const double &d) { if(!hasEnergies()) return; energies_->eslip_=d;}
        double  edashpot() const {return hasEnergies() ? energies_->edashpot_: 0.0;}
        void    edashpot(const double &d) { if(!hasEnergies()) return; energies_->edashpot_=d;}
        double  epbstrain() const {return hasEnergies() ? energies_->epbstrain_: 0.0;}
        void    epbstrain(const double &d) { if(!hasEnergies()) return; energies_->epbstrain_=d;}

        bool     hasPBond() const {return pbProps_ ? true:false;}
        int      pbState()   const {return hasPBond() ? pbProps_->pb_state_: 0;}
        void     pbState(int i) { if(!hasPBond()) return; pbProps_->pb_state_=i;}
        double   pbRmul() const {return (hasPBond() ? (pbProps_->pb_rmul_) : 0.0);}
        void     pbRmul(const double &d) { if(!hasPBond()) return; pbProps_->pb_rmul_=d;}
        double   pbKn() const {return (hasPBond() ? (pbProps_->pb_kn_) : 0.0);}
        void     pbKn(const double &d) { if(!hasPBond()) return; pbProps_->pb_kn_=d;}
        double   pbKs() const {return (hasPBond() ? (pbProps_->pb_ks_) : 0.0);}
        void     pbKs(const double &d) { if(!hasPBond()) return; pbProps_->pb_ks_=d;}
		double   pbKn1() const { return (hasPBond() ? (pbProps_->pb_kn1_) : 0.0); }
		void     pbKn1(const double &d) { if (!hasPBond()) return; pbProps_->pb_kn1_ = d; }
		double   pbKs1() const { return (hasPBond() ? (pbProps_->pb_ks1_) : 0.0); }
		void     pbKs1(const double &d) { if (!hasPBond()) return; pbProps_->pb_ks1_ = d; }
        double   pbMCF() const {return (hasPBond() ? (pbProps_->pb_mcf_) : 0.0);}
        void     pbMCF(const double &d) { if(!hasPBond()) return; pbProps_->pb_mcf_=d;}
        double   pbTen() const {return (hasPBond() ? (pbProps_->pb_ten_) : 0.0);}
        void     pbTen(const double &d) { if(!hasPBond()) return; pbProps_->pb_ten_=d;}
		double   pbCom() const { return (hasPBond() ? (pbProps_->pb_com_) : 0.0); }
		void     pbCom(const double &d) { if (!hasPBond()) return; pbProps_->pb_com_ = d; }
		//double   pbMb() const { return (hasPBond() ? (pbProps_->pb_mb_) : 0.0); }
		//void     pbMb(const double &d) { if (!hasPBond()) return; pbProps_->pb_mb_ = d; }
		double   pbTen1() const { return (hasPBond() ? (pbProps_->pb_ten1_) : 0.0); }
		void     pbTen1(const double &d) { if (!hasPBond()) return; pbProps_->pb_ten1_ = d; }
		double   pbCom1() const { return (hasPBond() ? (pbProps_->pb_com1_) : 0.0); }
		void     pbCom1(const double &d) { if (!hasPBond()) return; pbProps_->pb_com1_ = d; }
		//double   pbMb1() const { return (hasPBond() ? (pbProps_->pb_mb1_) : 0.0); }
		//void     pbMb1(const double &d) { if (!hasPBond()) return; pbProps_->pb_mb1_ = d; }
        double   pbCoh() const {return (hasPBond() ? (pbProps_->pb_coh_) : 0.0);}
        void     pbCoh(const double &d) { if(!hasPBond()) return; pbProps_->pb_coh_=d;}
		double   pbCoh1() const { return (hasPBond() ? (pbProps_->pb_coh1_) : 0.0); }
		void     pbCoh1(const double &d) { if (!hasPBond()) return; pbProps_->pb_coh1_ = d; }
        //double   pbFA() const {return (hasPBond() ? (pbProps_->pb_fa_) : 0.0);}
        //void     pbFA(const double &d) { if(!hasPBond()) return; pbProps_->pb_fa_=d;}
        DVect    pbF() const {return hasPBond() ? pbProps_->pb_F_: DVect(0.0);}
        void     pbF(const DVect &f) { if(!hasPBond()) return; pbProps_->pb_F_=f;}
        DAVect   pbM() const {return hasPBond() ? pbProps_->pb_M_: DAVect(0.0);}
        void     pbM(const DAVect &m) { if(!hasPBond()) return; pbProps_->pb_M_=m;}
		
		double   Mn() const { return (hasPBond() ? (pbProps_->mn_) : 0.0); }
		void     Mn(const double &d) { if (!hasPBond()) return; pbProps_->mn_ = d; }
		double   Ms() const { return (hasPBond() ? (pbProps_->ms_) : 0.0); }
		void     Ms(const double &d) { if (!hasPBond()) return; pbProps_->ms_ = d; }
		double   Mm() const { return (hasPBond() ? (pbProps_->mm_) : 0.0); }
		void     Mm(const double &d) { if (!hasPBond()) return; pbProps_->mm_ = d; }
		double   Yield() const { return (hasPBond() ? (pbProps_->Fyield_) : 0.0); }
		void     Yield(const double &d) { if (!hasPBond()) return; pbProps_->Fyield_ = d; }
		double   Ucn() const { return (hasPBond() ? (pbProps_->ucn_) : 0.0); }
		void     Ucn(const double &d) { if (!hasPBond()) return; pbProps_->ucn_ = d; }
		double   Ucn1() const { return (hasPBond() ? (pbProps_->ucn1_) : 0.0); }
		void     Ucn1(const double &d) { if (!hasPBond()) return; pbProps_->ucn1_ = d; }
		double   Ucs() const { return (hasPBond() ? (pbProps_->ucs_) : 0.0); }
		void     Ucs(const double &d) { if (!hasPBond()) return; pbProps_->ucs_ = d; }
		double   Ucb() const { return (hasPBond() ? (pbProps_->ucb_) : 0.0); }
		void     Ucb(const double &d) { if (!hasPBond()) return; pbProps_->ucb_ = d; }

		//double   Alpha() const { return (hasPBond() ? (pbProps_->alpha_) : 0.0); }
		//void     Alpha(const double& d) { if (!hasPBond()) return; pbProps_->alpha_ = d; }
		//double   Timing() const { return (hasPBond() ? (pbProps_->timing_) : 0.0); }
		//void     Timing(const double& d) { if (!hasPBond()) return; pbProps_->timing_ = d; }


		double   Dvar() const { return (hasPBond() ? (pbProps_->D_) : 0.0); }
		void     Dvar(const double &d) { if (!hasPBond()) return; pbProps_->D_ = d; }
		double   Fyi() const { return (hasPBond() ? (pbProps_->Fy_) : 0.0); }
		void     Fyi(const double &d) { if (!hasPBond()) return; pbProps_->Fy_ = d; }

		//the following properties are not in the keyword property list
		const double   & Dellam() const { return delta_lamda_; }
		void             Dellam(const double &d) { delta_lamda_ = d; }
		const double   & Prev_dvar() const { return prev_D_; }
		void             Prev_dvar(const double &d) { prev_D_ = d; }
		const double   & Prev_dellam() const { return prev_dellam_; }
		void             Prev_dellam(const double &d) { prev_dellam_ = d; }
		const double   & Normal_total_disp() const { return normal_total_disp_; }
		void             Normal_total_disp(const double &d) { normal_total_disp_ = d; }
		const double   & Shear_total_dispy() const { return shear_total_dispy_; }
		void             Shear_total_dispy(const double &d) { shear_total_dispy_ = d; }
		const double   & Shear_total_dispz() const { return shear_total_dispz_; }
		void             Shear_total_dispz(const double &d) { shear_total_dispz_ = d; }
		const double   & Shear_total_disp() const { return shear_total_disp_; }
		void             Shear_total_disp(const double &d) { shear_total_disp_ = d; }
		const double   & Bending_total_roy() const { return bending_total_roy_; }
		void             Bending_total_roy(const double &d) { bending_total_roy_ = d; }
		const double   & Bending_total_roz() const { return bending_total_roz_; }
		void             Bending_total_roz(const double &d) { bending_total_roz_ = d; }
		const double   & Bending_total_ro() const { return bending_total_ro_; }
		void             Bending_total_ro(const double &d) { bending_total_ro_ = d; }
		const double   & Normal_plastic_disp() const { return normal_plastic_disp_; }
		void             Normal_plastic_disp(const double &d) { normal_plastic_disp_ = d; }
		const double   & Shear_plastic_disp() const { return shear_plastic_disp_; }
		void             Shear_plastic_disp(const double &d) { shear_plastic_disp_ = d; }
		const double   & Bending_plastic_ro() const { return bending_plastic_ro_; }
		void             Bending_plastic_ro(const double &d) { bending_plastic_ro_ = d; }

		const double   & BendingMoy() const { return BendingMoy_; }
		void             BendingMoy(const double &d) { BendingMoy_ = d; }
		const double   & BendingMoz() const { return BendingMoz_; }
		void             BendingMoz(const double &d) { BendingMoz_ = d; }
		const double   & BendingMo() const { return BendingMo_; }
		void             BendingMo(const double &d) { BendingMo_ = d; }
		const double   & Vy() const { return  Vy_; }
		void             Vy(const double &d) { Vy_ = d; }
		const double   & Vz() const { return Vz_; }
		void             Vz(const double &d) { Vz_ = d; }


		const double   & Prev_npd() const { return prev_npd_; }
		void             Prev_npd(const double &d) { prev_npd_ = d; }
		const double   & Prev_spd() const { return prev_spd_; }
		void             Prev_spd(const double &d) { prev_spd_ = d; }
		const double   & Prev_bpr() const { return prev_bpr_; }
		void             Prev_bpr(const double &d) { prev_bpr_ = d; }

		const double   & Npd1() const { return npd1_; }
		void             Npd1(const double &d) {npd1_ = d; }
		const double   & Spd1() const { return spd1_; }
		void             Spd1(const double &d) { spd1_ = d; }
		const double   & Bpr1() const { return bpr1_; }
		void             Bpr1(const double &d) { bpr1_ = d; }

		const double   & ReFnorm() const { return force_return_norm_; }
		void             ReFnorm(const double &d) { force_return_norm_ = d; }
		const double   & ReFtany() const { return force_return_tany_; }
		void             ReFtany(const double &d) { force_return_tany_ = d; }
		const double   & ReFtanz() const { return force_return_tanz_; }
		void             ReFtanz(const double &d) { force_return_tanz_ = d; }
		const double   & ReMby() const { return moment_return_by_; }
		void             ReMby(const double &d) { moment_return_by_ = d; }
		const double   & ReMbz() const { return moment_return_bz_; }
		void             ReMbz(const double &d) { moment_return_bz_ = d; }

        DVect2   pbTransStiff() const {return hasPBond() ? pbProps_->pbTransStiff_: DVect2(0.0);}
        void     pbTransStiff(const DVect2 &f) { if(!hasPBond()) return; pbProps_->pbTransStiff_=f;}
        DAVect   pbAngStiff() const {return hasPBond() ? pbProps_->pbAngStiff_: DAVect(0.0);}
        void     pbAngStiff(const DAVect &m) { if(!hasPBond()) return; pbProps_->pbAngStiff_=m;}
		

        uint inheritanceField() const {return inheritanceField_;}
        void inheritanceField(uint i) {inheritanceField_ = i;}

        const DVect2 & effectiveTranslationalStiffness()  const          {return effectiveTranslationalStiffness_;}
        void           effectiveTranslationalStiffness(const DVect2 &v ) {effectiveTranslationalStiffness_=v;}

        /// Return the total force that the contact model holds.
        virtual DVect    getForce(const IContactMechanical *) const;

        /// Return the total moment on 1 that the contact model holds
        virtual DAVect   getMomentOn1(const IContactMechanical *) const;

        /// Return the total moment on 1 that the contact model holds
        virtual DAVect   getMomentOn2(const IContactMechanical *) const;

    private:
        static int index_;

        struct Energies {
            Energies() : estrain_(0.0), eslip_(0.0), edashpot_(0.0), epbstrain_(0.0) {}
            double estrain_;     // elastic energy stored in contact
            double eslip_;       // work dissipated by friction
            double edashpot_;    // work dissipated by dashpots
            double epbstrain_;   // parallel bond strain energy
        };

        struct dpProps {
            dpProps() : dp_nratio_(0.0), dp_sratio_(0.0), dp_mode_(0), dp_F_(DVect(0.0)) {}
            double dp_nratio_;     // normal viscous critical damping ratio
            double dp_sratio_;     // shear  viscous critical damping ratio
            int    dp_mode_;       // for viscous mode (0-4) 0 = dashpots, 1 = tensile limit, 2 = shear limit, 3 = limit both
            DVect  dp_F_;          // Force in the dashpots
        };

        struct pbProps {
			pbProps() : pb_state_(0), pb_rmul_(1.0), pb_kn_(0.0), pb_ks_(0.0), pb_kn1_(0.0), pb_ks1_(0.0),
				pb_mcf_(1.0), pb_ten_(0.0), pb_com_(0.0), pb_ten1_(0.0), pb_com1_(0.0), pb_coh_(0.0), pb_coh1_(0.0), pb_F_(DVect(0.0)), pb_M_(DAVect(0.0)),
				mn_(0.0), ms_(0.0), mm_(0.0), Fyield_(0.0), ucn_(0.0), ucn1_(0.0), ucs_(0.0), ucb_(0.0), /*alpha_(0.0), timing_(0.0),*/ D_(0.0), Fy_(0.0), pbTransStiff_(0.0), pbAngStiff_(0.0) {}
            // parallel bond
            int     pb_state_;                            // Bond mode
            double  pb_rmul_;                             // Radius multiplier
            double  pb_kn_;                               // initial normal stiffness
            double  pb_ks_;                               // initial shear stiffness
			double  pb_kn1_;                              // damaged normal stiffness
			double  pb_ks1_;                              // damaged shear stiffness
            double  pb_mcf_;                              // Moment contribution factor
            double  pb_ten_;                              // initial tensile strength
			double  pb_com_;                              // initial compressive strength
			//double  pb_mb_;                               // initial bending moment strength
			double  pb_ten1_;                             // damaged tensile strength
			double  pb_com1_;                             // damaged compressive strength
			//double  pb_mb1_;                              // damaged bending moment strength
			double  pb_coh_;                              // initial cohesion
			double  pb_coh1_;                             // damaged cohesion
            //double  pb_fa_;                               // friction angle
			DVect   pb_F_;                                // Force in parallel bond
            DAVect  pb_M_;                                // moment in parallel bond
			double  mn_;                                  // normal plastic return direction i.e. the normal slope in the yield surface
			double  ms_;                                  // shear plastic return direction i.e. the shear slope in the yield surface
			double  mm_;                                  // bending moment return direction
			double  Fyield_;                              // value of the yield function
			double  ucn_;                                 // max value of normal plastic disp 
			double  ucn1_;                                // softening for compression
			double  ucs_;                                 // max value of shear plastic disp 
			double  ucb_;                                 // max value of bending angle 
			//double  alpha_;                                 
			//double  timing_;
			double  D_;                                   // damage variable
			double  Fy_;
            DVect2  pbTransStiff_;    // (Normal,Shear) Translational stiffness of the parallel bond
            DAVect  pbAngStiff_;      // (Normal,Shear) Rotational stiffness of the parallel bond
        };



        bool   updateKn(const IContactMechanical *con);
        bool   updateKs(const IContactMechanical *con);
        bool   updateFric(const IContactMechanical *con);
        double pbStrainEnergy() const;                     // Compute  bond strain energy

        void   updateEffectiveStiffness(ContactModelMechanicalState *state);

        DVect3 pbData(const IContactMechanical *con) const; // Bond area and inertia
        DVect2 pbSMax(const IContactMechanical *con) const; // Maximum stress (tensile,shear) at bond periphery
        //double pbShearStrength(const double &pbArea) const;      // Bond shear strength
        void   setDampCoefficients(const double &mass,double *vcn,double *vcs);

        // inheritance fields
        quint32 inheritanceField_;

        // linear model
        double      kn_;        // normal stiffness
        double      ks_;        // shear stiffness
        double      fric_;      // Coulomb friction coefficient
        DVect       lin_F_;     // Force carried in the linear model
        bool        lin_S_;     // the current sliding state
        uint        lin_mode_;  // specifies absolute (0) or incremental (1) behavior for the the linear part
        double      rgap_;      // reference gap for the linear part

        dpProps *   dpProps_;     // The viscous properties
        pbProps *   pbProps_;     // The parallel bond properties
        double      userArea_;    // User specified area
		
		double  delta_lamda_;                         // plastic multiplier
		double  prev_D_;                              // damage variable in the previous step
		double  prev_dellam_;                         // plastic multiplier in the previous step
		double  normal_total_disp_;                   // normal total disp
		double  shear_total_dispy_;
		double  shear_total_dispz_;
		double  shear_total_disp_;                    // shear total disp
		double  bending_total_roy_;
		double  bending_total_roz_;
		double  bending_total_ro_;                    // bending total ro
		double  normal_plastic_disp_;                 // normal plastic disp
		double  shear_plastic_disp_;                  // shear plastic disp
		double  bending_plastic_ro_;                  // bending plastic ro
		double  BendingMoy_;
		double  BendingMoz_;
		double  BendingMo_;
		double  Vy_;
		double  Vz_;
		double  prev_npd_;                            // the normal plastic disp at timestep k-2
		double  prev_spd_;
		double  prev_bpr_;
		double  npd1_;                            // the normal plastic disp at timestep k-2
		double  spd1_;
		double  bpr1_;
		double  force_return_norm_;                   // normal force return
		double  force_return_tany_;                   // shear force return in y direction
		double  force_return_tanz_;                   // shear force return in z direction
		double  moment_return_by_;                    // bening moment return in y direction
		double  moment_return_bz_;                    // bening moment return in z direction

        Energies *  energies_;    // energies

        DVect2  effectiveTranslationalStiffness_;
		
    };
} // namespace itascaxd
// EoF



