FC = gfortran

FFLAGS = -O3 -ffree-form
# debugging
#FFLAGS = -g -ffree-form

HYDRO_FFLAGS = -ffree-line-length-none -I/usr/local/include


# modules
BEHAVIOR_MOD = behavior_module.o

BOUNDARY_MOD = boundary_module.o

CONVERT_MOD = conversion_module.o

GRIDCELL_MOD = gridcell_module.o

HTURB_MOD = hor_turb_module.o

HYDRO_MOD = hydrodynamic_module.o

INT_MOD = interpolation_module.o

PARAM_MOD = parameter_module_wb.o

NORM_MOD = norm_module.o

PIP_MOD = point_in_polygon_module.o

RANDOM_MOD = random_module.o

SETTLEMENT_MOD = settlement_module.o

TENSION_MOD = tension_module.o

VTURB_MOD = ver_turb_module.o


OBJ = $(BEHAVIOR_MOD) $(BOUNDARY_MOD) $(CONVERT_MOD) $(GRIDCELL_MOD) $(HTURB_MOD) $(HYDRO_MOD) $(INT_MOD) $(PARAM_MOD) $(NORM_MOD) $(PIP_MOD) $(RANDOM_MOD) $(SETTLEMENT_MOD) $(TENSION_MOD) $(VTURB_MOD)

LIB = -L/usr/local/lib -lnetcdf -lhdf5_hl -lhdf5 /usr/local/lib/libnetcdff.a /usr/local/lib/libnetcdf.a

.SUFFIXES : .f90

LTRANS.exe: $(OBJ) LTRANS.f90
	$(FC) LTRANS.f90 -o $@ $(FFLAGS) $(OBJ) $(LIB)

$(BEHAVIOR_MOD): $(PARAM_MOD) $(SETTLEMENT_MOD) $(HYDRO_MOD) $(NORM_MOD) $(RANDOM_MOD)

$(BOUNDARY_MOD): $(PARAM_MOD) $(HYDRO_MOD) $(PIP_MOD) $(CONVERT_MOD)

$(CONVERT_MOD): $(PARAM_MOD)

$(HTURB_MOD): $(PARAM_MOD) $(NORM_MOD)

$(HYDRO_MOD): $(GRIDCELL_MOD) $(PARAM_MOD) $(TENSION_MOD) $(INT_MOD) $(CONVERT_MOD) hydrodynamic_module.f90
	$(FC) $(FFLAGS) $(HYDRO_FFLAGS) -c hydrodynamic_module.f90

$(PARAM_MOD): GRID_wb.inc LTRANS_wb.inc

$(NORM_MOD): $(RANDOM_MOD) $(PARAM_MOD)

$(SETTLEMENT_MOD): $(PARAM_MOD) $(CONVERT_MOD) $(HYDRO_MOD) $(GRIDCELL_MOD) $(PIP_MOD)

$(VTURB_MOD): $(PARAM_MOD) $(HYDRO_MOD) $(INT_MOD) $(NORM_MOD) $(TENSION_MOD)

.f90.o:
	$(FC) $(FFLAGS) -c $<

wb:
	./LTRANS_compile_wb.sh

clean_output:
	rm -rf para*.csv endfile.csv

clean: clean_output
	rm -rf LTRANS.exe *.o *.mod
