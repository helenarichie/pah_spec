# intel compiler on ourania5
#FC=ifort -c -O2
#FL=ifort

# debug:
#FC=ifort -c -C
#FL=ifort -C

# gfortran
FC=gfortran -c -O2 -C
FL=gfortran

VPATH = ./src

OBJDIR = ./obj
EXEDIR = ./exe

TDISTOBJS = vsg_stat_therm.o agli_fzero.o agli_zbrent_d_v2.o \
	bb_radiation_field.o bicg.o \
	bruzual_charlot_rf.o cooling_rate.o \
	dbgdl_cooling_afi_v2.o debye.o dens_states.o dqg32.o \
        ext_hd20.o fe_degeneracy.o fe_density_states_v2.o \
	fe_enthalpy_state.o fe_modes.o fe_spec_heat_v1.o \
	fe_temp.o glu.o gra_pah_degeneracy.o \
	gra_pah_density_states.o gra_pah_enthalpy_state.o \
	heating_afi_v2.o heating_rate_v2.o hpsort_d.o \
	mmp_radiation_field.o \
	pah_modes.o pah_size_atoms.o pah_spec_heat.o \
	pah_temp.o parab_interpl.o planck.o qcomp_pack.o \
	readqlib.o sil_degeneracy.o sil_density_states.o \
	sil_enthalpy_state_v2.o sil_gra_pah_cs_v6.o \
	sil_gra_pah_temp_v11.o sil_modes.o sil_spec_heat.o \
	sil_temp.o silprof_gc.o simpson_dbgdl_cooling_rate_v2.o \
	simpson_heating_rate_v2.o simpson_intrabin_heating_rate.o \
	simpson_intrabin_therm_cooling_rate_v2.o \
	simpson_stat_cooling_rate_v2.o simpson_therm_cooling_rate_v2.o \
	spline.o splint.o starburst_radiation_field.o \
	stat_cooling_afi_v2.o therm_cooling_afi_v2.o \
	transition_matrix_v4.o trap_integ.o vsg_td_emission_v12.o \
	zbrent_d.o

vsg_stat_therm: $(addprefix $(OBJDIR)/, $(TDISTOBJS))
	$(FL) $^ $(LDLIBS) -o $(addprefix $(EXEDIR)/, $@)

$(OBJDIR)/%.o: %.f
	${FC} $< -o $@

clean:
	$(RM) $(OBJDIR)/*.o
	$(RM) .depend
