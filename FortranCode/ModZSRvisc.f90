module ModZSRvisc

use Mod_kind_param, only : wp
use ModSystemProp, only : nneutral, calcviscosity, nindcomp
use ModAIOMFACvar, only : activity, ln_etamix, ln_etamixZSR, XrespSalt, wtf, T_K, ln_etaZSR_org, ln_etaZSR_inorg
    
implicit none
!public variables:
real(wp),dimension(2),public :: ZSRviscmix
!module-private variables:
real(wp),dimension(:),allocatable,private :: WTFsub, WTFsub_old
real(wp),private :: aw_fullsys, drymass

!public procedures:
public :: ZSRviscosity
private     !default
    
!make all variables of this module threadprivate for use in parallel execution with openMP:
!$OMP THREADPRIVATE( WTFsub, WTFsub_old, aw_fullsys, drymass, ZSRviscmix )

!==============!
contains 
!==============!

    
    !-------------------------------------------------------------------------------------
    subroutine ZSRviscosity(wtf_, TK_)

    use ModCalcActCoeff, only : AIOMFAC_calc
    use Mod_MINPACK, only : hybrd1

    implicit none
    !interface vars
    real(wp),dimension(:),intent(in) :: wtf_
    real(wp),intent(in) :: TK_
    !local vars
    integer,parameter :: n = 1
    integer :: ii, info, iflag
    real(wp),parameter :: deps = epsilon(1.0_wp), sqrt_deps = sqrt(deps), sqrt_tiny = sqrt(tiny(1.0_wp))
    real(wp) :: OIR, f1, f2, sumwtforg, sumwtfinorg, tol
    real(wp),dimension(n) :: fvec_inorg, fvec_org, m1
    real(wp),dimension(nindcomp) :: wtf_fullsys, WTFsub_inorg, WTFsub_org, XrespSalt_fullsys 
    logical,parameter :: onlyDeltaVisc = .true.
    logical :: calcviscosity_saved
    !........................................
        
    calcviscosity_saved = calcviscosity !AZ added this; 2021-11-24
    calcviscosity = .false.
    
    call AIOMFAC_calc(wtf_, TK_)        !calculate full system water activity
    aw_fullsys = activity(1)            !save full system water activity
    if (aw_fullsys < sqrt_tiny) then
        write(*,*) "aw_fullsys ", aw_fullsys
        write(*,*) "wtf_ ", wtf_
        aw_fullsys = sqrt_tiny
    endif
    if (aw_fullsys < 5.0E-2_wp) then
        tol = sqrt_deps
    else
        tol = 1.0E2_wp*sqrt_deps
    endif    
    wtf_fullsys = wtf
    XrespSalt_fullsys = XrespSalt
    sumwtforg = sum(wtf_fullsys(2:nneutral))
    sumwtfinorg = sum(wtf_fullsys(nneutral+1:))
    
    !if ZSRviscosity is called for a mixture without organic or without inorganic
    if (sumwtfinorg > 0.0_wp) then
        OIR = sumwtforg/sumwtfinorg
    else
        OIR = 1.0_wp
    endif
    
    !inorganic subsystem
    allocate( WTFsub(nindcomp), WTFsub_old(nindcomp) )
    WTFsub = wtf_fullsys
    WTFsub(2:nneutral) = 0.0_wp              !zero out the organics
    WTFsub = WTFsub/sum(WTFsub)
    WTFsub_old = WTFsub                     !this will have the old version of the subsystem...
    drymass = sum(WTFsub_old(2:))
    !use the hybrd1 method from Mod_MINPACK to solve an equation (or as system of eqs.) numerically; here n = 1
    do ii = 1,5
        !generate a good initial guess for electrolyte-associated water mass:
        m1(1) = wtf_fullsys(1)**(1.0_wp/ii) * sumwtfinorg/(sumwtforg + sumwtfinorg)
        !---
        call hybrd1 ( fcn1ZSRvisc, n, m1, fvec_inorg, tol, info )
        !---
        if (info == 1) then     !this should usually work
            exit
        else
            iflag = info
        endif
    enddo
    
    !re-compute viscosity of the inorganic subsystem with the determined water amount:
    calcviscosity = .true.
    call fcn1ZSRvisc( n, m1, fvec_inorg, iflag )
    WTFsub_inorg = WTFsub
    ln_etaZSR_inorg = ln_etamix
    
    !organic subsystem
    calcviscosity = .false.
    WTFsub = wtf_fullsys
    WTFsub(nneutral+1:) = 0.0_wp             !zero out the electrolytes
    WTFsub = WTFsub/sum(WTFsub)
    WTFsub_old = WTFsub                     !this will have the old version of the subsystem... 
    drymass = sum(WTFsub_old(2:))
    do ii = 1,5
        !generate a good initial guess for organics-associated water mass:
        m1(1) = wtf_fullsys(1)**(1.0_wp/ii) * sumwtforg/(sumwtforg + sumwtfinorg)
        !---
        call hybrd1( fcn1ZSRvisc, n, m1, fvec_org, tol, info )
        !---
        if (info == 1) then
            exit
        else
            iflag = info
        endif
    enddo
    
    !re-compute viscosity of the organic subsystem with the determined water amount:
    calcviscosity = .true.
    call fcn1ZSRvisc( n, m1, fvec_org, iflag )
    WTFsub_org = WTFsub 
    ln_etaZSR_org = ln_etamix
    
    !if ZSRviscosity is called for a mixture without organic or without inorganic
    if (sumwtforg < deps) then
        f1 = 0.0_wp
        f2 = 1.0_wp
        ZSRviscmix(1) = f1
        ZSRviscmix(2) = f2
    elseif (sumwtfinorg < deps) then
        f1 = 1.0_wp
        f2 = 0.0_wp
        ZSRviscmix(1) = f1
        ZSRviscmix(2) = f2
    else
        f1 = ( OIR * sum(WTFsub_inorg(nneutral+1:)) )/( sum(WTFsub_org(2:nneutral)) + OIR * sum(WTFsub_inorg(nneutral+1:)) )
        f2 = 1.0_wp - f1
        ZSRviscmix(1) = f1
        ZSRviscmix(2) = f2
    endif
    !*@#! ln_etamixZSR = exp(f1*log(etaZSR_org) + f2*log(etaZSR_inorg))
    ln_etamixZSR = f1*ln_etaZSR_org + f2*ln_etaZSR_inorg
    ln_etamix = ln_etamixZSR
    
    calcviscosity = calcviscosity_saved     ![changed by AZ; 2021-11-24 from calcviscosity = .false.] 
    
    deallocate(WTFsub, WTFsub_old)

    wtf = wtf_fullsys

    end subroutine ZSRviscosity
    !-------------------------------------------------------------------------------------

    
    !-------------------------------------------------------------------------------------
    subroutine fcn1ZSRvisc( n, watermass, fvec, iflag )

    use ModCalcActCoeff, only : AIOMFAC_calc

    implicit none

    !interface vars:
    integer,intent(in) :: n
    integer,intent(out) :: iflag
    real(wp),dimension(n),intent(inout) :: watermass
    real(wp),dimension(n),intent(out) :: fvec
    !local vars:
    real(wp),parameter :: deps = epsilon(1.0_wp), tinym1 = 1.0E3_wp*deps
    real(wp) :: awratio
    !.......................................
    iflag = 0
    
    !Constrain watermass value to be positive:
    watermass(1) = max(tinym1, watermass(1))
            
    WTFsub(1) = watermass(1) / (watermass(1) + drymass)
    WTFsub(2:) = WTFsub_old(2:) / (watermass(1) + drymass)      !scale the dry mass to the new water content

    call AIOMFAC_calc(WTFsub, T_K)
    
    !compute relative deviation as current water activity from the subsystem over the saved aw from full system:
    awratio = activity(1)/aw_fullsys
    if (awratio < deps) then        !deal with numerical precision limitations
        fvec(1) = log(awratio)      !this fvec is not the relative deviation; however, this expression should allow the solver to recover such that awratio becomes > deps again; 
    else                            !the usual case
        fvec(1) = awratio - 1.0_wp
    endif

    end subroutine fcn1ZSRvisc
    !-------------------------------------------------------------------------------------

end module ModZSRvisc