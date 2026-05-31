/*
    page_registry.hpp
    Registry of all main-menu pages.

    HOW TO ADD A NEW PAGE
    ─────────────────────
    1. Create  src/MAIN-MENU/pages/my_page.hpp
       • Define  void draw(AppState& state)  inside a unique namespace, e.g. pages::mypage
    2. Add a new entry to the Tab enum in app_state.hpp
    3. #include your header below, then append one PageDef to kPages.

    That's it, the sidebar and content area are driven by this list automatically.
*/
#pragma once

#include "app_state.hpp"

// ── page headers ──────────────────────────────────────────────────────────────
#include "pages/overview.hpp"
#include "pages/simulation.hpp"
#include "pages/data_analysis.hpp"
#include "pages/object_library.hpp"
#include "pages/export.hpp"
#include "pages/settings.hpp"
#include "pages/breakdown.hpp"
// ─────────────────────────────────────────────────────────────────────────────

#include <span>

namespace ui {

    struct PageDef {
        const char* label;    // sidebar label
        const char* section;  // sidebar group heading (nullptr = same section as previous)
        Tab         id;
        void (*draw)(AppState&);
    };

    // ── THE PAGE LIST ─────────────────────────────────────────────────────────
    // Sections are printed as headings the first time they appear.
    inline constexpr PageDef kPages[] = {
        { "Overview",       "WORKSPACE", Tab::Overview,      pages::overview::draw      },
        { "Simulation",     nullptr,     Tab::Simulation,    pages::simulation::draw    },
        { "Data Analysis",  nullptr,     Tab::DataAnalysis,  pages::data_analysis::draw },
        { "Object Library", nullptr,     Tab::ObjectLibrary, pages::object_library::draw},
        { "Export",         "TOOLS",     Tab::Export,        pages::export_tab::draw    },
        { "Settings",       nullptr,     Tab::Settings,      pages::settings::draw      },
        { "Breakdown",      nullptr,     Tab::Breakdown,     pages::breakdown::draw     },
    };
    // ─────────────────────────────────────────────────────────────────────────

    inline std::span<const PageDef> getPages() {
        return std::span<const PageDef>(kPages);
    }

} // namespace ui
