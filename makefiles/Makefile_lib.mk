-include project_config.mk
-include dependencies.mk
-include $(ROOT_DIR)/makefiles/recursive_deps.mk
-include $(ROOT_DIR)/makefiles/common.mk

PROJECT_NAME=$(subst $(ROOT_DIR)/src/,,$(shell pwd))
SRCS=$(wildcard *.cpp)
OBJS=$(patsubst %.cpp,$(TEMPDIR)/$(PROJECT_NAME)/%.o,$(SRCS))
DEBUG_OBJS=$(patsubst %.cpp,$(TEMPDIR)/debug/$(PROJECT_NAME)/%.o,$(SRCS))
PROF_OBJS=$(patsubst %.cpp,$(TEMPDIR)/prof/$(PROJECT_NAME)/%.o,$(SRCS))

TARGET=$(BINDIR)/$(PROJECT_NAME).a
DEBUG_TARGET=$(BINDIR)/debug/$(PROJECT_NAME).a
PROF_TARGET=$(BINDIR)/prof/$(PROJECT_NAME).a

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

$(TARGET): $(OBJS)
	@mkdir -p $(dir $@)
	$(call log_ar,$(PROJECT_NAME),$(subst $(ROOT_DIR)/,,$@))
	@$(AR) rcs $@ $(OBJS)
	@echo ""

$(DEBUG_TARGET): $(DEBUG_OBJS)
	@mkdir -p $(dir $@)
	$(call log_ar,$(PROJECT_NAME),$(subst $(ROOT_DIR)/,,$@))
	@$(AR) rcs $@ $(DEBUG_OBJS)
	@echo ""

$(PROF_TARGET): $(PROF_OBJS)
	@mkdir -p $(dir $@)
	$(call log_ar,$(PROJECT_NAME),$(subst $(ROOT_DIR)/,,$@))
	@$(AR) rcs $@ $(PROF_OBJS)
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

-include $(OBJS_DEPEND)
-include $(DEBUG_OBJS_DEPEND)
-include $(PROF_OBJS_DEPEND)
