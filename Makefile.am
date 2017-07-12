AUTOMAKE_OPTIONS = check-news
ACLOCAL_AMFLAGS = -I local_m4

MNI_MODELS := mni-models

# mess with the data install directory
pkgdatadir = $(prefix)/share/$(MNI_MODELS)

MODEL_BASE := mni_icbm152
T1 := t1
T2 := t2
PD := pd
CLS := cls
CSF := csf
GM := gm
WM := wm
PRIM_MODAL := $(T1)
SYM_MODEL := tal_nlin_sym_09a
ASYM_MODEL := tal_nlin_asym_09a
PRIMARY_MODELS := $(MODEL_BASE)_$(PRIM_MODAL)_$(SYM_MODEL) \
		  $(MODEL_BASE)_$(PRIM_MODAL)_$(ASYM_MODEL)
SECONDARY_MODELS := $(MODEL_BASE)_$(T2)_$(SYM_MODEL) \
		    $(MODEL_BASE)_$(T2)_$(ASYM_MODEL) \
		    $(MODEL_BASE)_$(PD)_$(SYM_MODEL) \
		    $(MODEL_BASE)_$(PD)_$(ASYM_MODEL)

# datestamp
DATE_STAMP := `date +%D-%T`

# no-install package data

noinst_DATA := \
	$(PRIMARY_MODELS:=.mnc.gz) \
	$(PRIMARY_MODELS:=_mask.mnc.gz) \
	$(PRIMARY_MODELS:=_eye_mask.mnc.gz) \
	$(PRIMARY_MODELS:=_mask.mnc.gz) \
	$(PRIMARY_MODELS:=_mask.obj.gz) \
	$(PRIMARY_MODELS:=_face_mask.mnc.gz) \
	$(PRIMARY_MODELS:=_headmask.mnc.gz) \
	$(PRIMARY_MODELS:=_vent.mnc.gz) \
	$(PRIMARY_MODELS:=_vent.obj.gz) \
	$(PRIMARY_MODELS:=_outline.mnc.gz) \
	$(SECONDARY_MODELS:=.mnc.gz) \
	$(MODEL_BASE)_$(CLS)_$(SYM_MODEL).mnc.gz \
	$(MODEL_BASE)_$(CLS)_$(ASYM_MODEL).mnc.gz \
	$(MODEL_BASE)_$(CSF)_$(SYM_MODEL).mnc.gz \
	$(MODEL_BASE)_$(CSF)_$(ASYM_MODEL).mnc.gz \
	$(MODEL_BASE)_$(GM)_$(SYM_MODEL).mnc.gz \
	$(MODEL_BASE)_$(GM)_$(ASYM_MODEL).mnc.gz \
	$(MODEL_BASE)_$(WM)_$(SYM_MODEL).mnc.gz \
	$(MODEL_BASE)_$(WM)_$(ASYM_MODEL).mnc.gz \
	Cerebellum_Ventricles_SubCortical_Mask-2009s.mnc.gz

# ANIMAL files, should not really be here.
# mni_icbm152_lob_tal_nlin_sym_09a.mnc.gz

# installed data
pkgdata_DATA := \
	$(PRIMARY_MODELS:=.mnc) \
	$(PRIMARY_MODELS:=_mask.mnc) \
	$(PRIMARY_MODELS:=_mask.obj) \
	$(PRIMARY_MODELS:=_eye_mask.mnc) \
	$(PRIMARY_MODELS:=_face_mask.mnc) \
	$(PRIMARY_MODELS:=_headmask.mnc) \
	$(PRIMARY_MODELS:=_vent.mnc) \
	$(PRIMARY_MODELS:=_vent.obj) \
	$(PRIMARY_MODELS:=_outline.mnc) \
	$(SECONDARY_MODELS:=.mnc) \
	$(MODEL_BASE)_$(CLS)_$(SYM_MODEL).mnc \
	$(MODEL_BASE)_$(CLS)_$(ASYM_MODEL).mnc \
	$(MODEL_BASE)_$(CSF)_$(SYM_MODEL).mnc \
	$(MODEL_BASE)_$(CSF)_$(ASYM_MODEL).mnc \
	$(MODEL_BASE)_$(GM)_$(SYM_MODEL).mnc \
	$(MODEL_BASE)_$(GM)_$(ASYM_MODEL).mnc \
	$(MODEL_BASE)_$(WM)_$(SYM_MODEL).mnc \
	$(MODEL_BASE)_$(WM)_$(ASYM_MODEL).mnc \
	Cerebellum_Ventricles_SubCortical_Mask-2009s.mnc

EXTRA_DIST = $(noinst_DATA) $(m4_files) \
	epm-header.in

m4_files = local_m4/ac_prog_gzip.m4

# be tidy
CLEANFILES = $(pkgdata_DATA)

# gunzip rule
%: %.gz
	$(gzip) -d -c $< > $@

all: $(noinst_DATA:.gz=) $(DATE)
	for i in $(pkgdata_DATA); do \
		if [ $${i: -4} == ".mnc" ]; then \
		echo "Adding mni-model header to $$i"; \
		minc_modify_header \
			-sinsert mni-model:package=$(PACKAGE) \
			-sinsert mni-model:version=$(VERSION) \
			-sinsert mni-model:compile-date=$(DATE_STAMP) \
			-sinsert mni-model:compile-user=$$USER \
			$$i; \
		fi; \
	done
