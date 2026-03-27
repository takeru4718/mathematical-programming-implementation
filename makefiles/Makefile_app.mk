-include project_config.mk
-include dependencies.mk
-include $(ROOT_DIR)/makefiles/recursive_deps.mk
-include $(ROOT_DIR)/makefiles/common.mk

PROJECT_NAME=$(subst $(ROOT_DIR)/src/,,$(shell pwd))

SRCS=$(wildcard *.cpp)
OBJS=$(patsubst %.cpp,$(TEMPDIR)/$(PROJECT_NAME)/%.o,$(SRCS))
DEBUG_OBJS=$(patsubst %.cpp,$(TEMPDIR)/debug/$(PROJECT_NAME)/%.o,$(SRCS))
PROF_OBJS=$(patsubst %.cpp,$(TEMPDIR)/prof/$(PROJECT_NAME)/%.o,$(SRCS))

TARGET=$(BINDIR)/$(PROJECT_NAME)
DEBUG_TARGET=$(BINDIR)/debug/$(PROJECT_NAME)
PROF_TARGET=$(BINDIR)/prof/$(PROJECT_NAME)
DEPEND_DATA=$(patsubst $(ROOT_DIR)/data/$(PROJECT_NAME)/%,$(ROOT_DIR)/debug/$(PROJECT_NAME)/%,$(wildcard $(ROOT_DIR)/data/$(PROJECT_NAME)/*))

OBJS_DEPEND=$(OBJS:.o=.d)
DEBUG_OBJS_DEPEND=$(DEBUG_OBJS:.o=.d)
PROF_OBJS_DEPEND=$(PROF_OBJS:.o=.d)

LIB_OPTS=-L$(BINDIR) $(addprefix -l, $(DEPEND_LIBS))
DEBUG_LIB_OPTS=-L$(BINDIR)/debug $(addprefix -l, $(DEPEND_LIBS))
PROF_LIB_OPTS=-L$(BINDIR)/prof $(addprefix -l, $(DEPEND_LIBS))
INCLUDE_DIR_OPTS=$(addprefix -I$(ROOT_DIR)/src/lib, $(DEPEND_LIBS))

.PHONY: build
build: $(TARGET)
	@:
	
.PHONY: debug-build
debug-build: $(DEBUG_TARGET)
	@:

.PHONY: prof-build
prof-build: $(PROF_TARGET)
	@:

.PHONY: run
run: build $(DEPEND_DATA)
	@mkdir -p $(ROOT_DIR)/debug/$(PROJECT_NAME)
	$(call log_exec,$(PROJECT_NAME),$(subst $(ROOT_DIR)/,,$(TARGET)))
	@(cd $(ROOT_DIR)/debug/$(PROJECT_NAME) && $(TARGET) $(ARGS))

.PHONY: prof-run
prof-run: prof-build $(DEPEND_DATA)
	@mkdir -p $(ROOT_DIR)/debug/$(PROJECT_NAME)
	$(call log_exec,$(PROJECT_NAME),$(subst $(ROOT_DIR)/,,$(PROF_TARGET)))
	@(cd $(ROOT_DIR)/debug/$(PROJECT_NAME) && $(PROF_TARGET) $(ARGS))
	@gprof $(PROF_TARGET) $(ROOT_DIR)/debug/$(PROJECT_NAME)/gmon.out > $(ROOT_DIR)/debug/$(PROJECT_NAME)/gprof.txt
	@echo "Profile data saved to $(ROOT_DIR)/debug/$(PROJECT_NAME)/gprof.txt"

$(TARGET): $(OBJS) $(patsubst %, $(BINDIR)/lib%.a, $(DEPEND_LIBS))
	@mkdir -p $(dir $@)
	$(call log_ld,$(PROJECT_NAME),$(subst $(ROOT_DIR)/,,$@))
	@$(CXX) -o $@ $(OBJS) $(CXXFLAGS) $(LIB_OPTS)
	@echo ""

$(DEBUG_TARGET): $(DEBUG_OBJS) $(patsubst %, $(BINDIR)/debug/lib%.a, $(DEPEND_LIBS))
	@mkdir -p $(dir $@)
	$(call log_ld,$(PROJECT_NAME),$(subst $(ROOT_DIR)/,,$@))
	@$(CXX) -o $@ $(DEBUG_OBJS) $(CXXFLAGS) -O0 -g $(DEBUG_LIB_OPTS)
	@echo ""

$(PROF_TARGET): $(PROF_OBJS) $(patsubst %, $(BINDIR)/prof/lib%.a, $(DEPEND_LIBS))
	@mkdir -p $(dir $@)
	$(call log_ld,$(PROJECT_NAME),$(subst $(ROOT_DIR)/,,$@))
	@$(CXX) -o $@ $(PROF_OBJS) $(CXXFLAGS) -O0 -pg $(PROF_LIB_OPTS)
	@echo ""

$(OBJS): $(TEMPDIR)/$(PROJECT_NAME)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(call log_cxx,$(PROJECT_NAME),$<)
	@$(CXX) -c $< -o $@ $(CXXFLAGS) $(INCLUDE_DIR_OPTS) -MMD -MP

$(DEBUG_OBJS): $(TEMPDIR)/debug/$(PROJECT_NAME)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(call log_cxx,$(PROJECT_NAME),$<)
	@$(CXX) -c $< -o $@ $(CXXFLAGS) $(INCLUDE_DIR_OPTS) -O0 -g -MMD -MP

$(PROF_OBJS): $(TEMPDIR)/prof/$(PROJECT_NAME)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(call log_cxx,$(PROJECT_NAME),$<)
	@$(CXX) -c $< -o $@ $(CXXFLAGS) $(INCLUDE_DIR_OPTS) -O0 -pg -MMD -MP

$(patsubst %, $(BINDIR)/lib%.a, $(DEPEND_LIBS)): $(BINDIR)/%.a: FORCE
	@$(MAKE) -C $(ROOT_DIR) build-$*

$(patsubst %, $(BINDIR)/debug/lib%.a, $(DEPEND_LIBS)): $(BINDIR)/debug/%.a: FORCE
	@$(MAKE) -C $(ROOT_DIR) debug-build-$*

$(patsubst %, $(BINDIR)/prof/lib%.a, $(DEPEND_LIBS)): $(BINDIR)/prof/%.a: FORCE
	@$(MAKE) -C $(ROOT_DIR) prof-build-$*

$(DEPEND_DATA): $(ROOT_DIR)/debug/$(PROJECT_NAME)/%: $(ROOT_DIR)/data/$(PROJECT_NAME)/%
	@mkdir -p $(dir $@)
	@cp $< $@

.PHONY: clean
clean:
	$(call log_rm,$(PROJECT_NAME),$(subst $(ROOT_DIR)/,,$(TARGET)))
	@rm -f $(TARGET)

	$(call log_rm,$(PROJECT_NAME),$(subst $(ROOT_DIR)/,,$(DEBUG_TARGET)))
	@rm -f $(DEBUG_TARGET)

	$(call log_rm,$(PROJECT_NAME),$(subst $(ROOT_DIR)/,,$(PROF_TARGET)))
	@rm -f $(PROF_TARGET)

	$(call log_rm,$(PROJECT_NAME),$(subst $(ROOT_DIR),,$(TEMPDIR))/$(PROJECT_NAME))
	@rm -rf $(TEMPDIR)/$(PROJECT_NAME)

	$(call log_rm,$(PROJECT_NAME),$(subst $(ROOT_DIR),,$(TEMPDIR))/debug/$(PROJECT_NAME))
	@rm -rf $(TEMPDIR)/debug/$(PROJECT_NAME)

	$(call log_rm,$(PROJECT_NAME),$(subst $(ROOT_DIR),,$(TEMPDIR))/prof/$(PROJECT_NAME))
	@rm -rf $(TEMPDIR)/prof/$(PROJECT_NAME)

.PHONY: debug-clean
debug-clean:
	$(call log_rm,$(PROJECT_NAME),debug/$(PROJECT_NAME))
	@rm -rf $(ROOT_DIR)/debug/$(PROJECT_NAME)

.PHONY: FORCE
FORCE: ;

-include $(OBJS_DEPEND)
-include $(DEBUG_OBJS_DEPEND)
-include $(PROF_OBJS_DEPEND)
