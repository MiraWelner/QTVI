#include "noise_marking_gui.h"
#include "noise_marking_guiPlugin.h"

#include <QtCore/QtPlugin>

noise_marking_guiPlugin::noise_marking_guiPlugin(QObject *parent)
    : QObject(parent)
{
    initialized = false;
}

void noise_marking_guiPlugin::initialize(QDesignerFormEditorInterface * /*core*/)
{
    if (initialized)
        return;

    initialized = true;
}

bool noise_marking_guiPlugin::isInitialized() const
{
    return initialized;
}

QWidget *noise_marking_guiPlugin::createWidget(QWidget *parent)
{
    return new noise_marking_gui(parent);
}

QString noise_marking_guiPlugin::name() const
{
    return "noise_marking_gui";
}

QString noise_marking_guiPlugin::group() const
{
    return "My Plugins";
}

QIcon noise_marking_guiPlugin::icon() const
{
    return QIcon();
}

QString noise_marking_guiPlugin::toolTip() const
{
    return QString();
}

QString noise_marking_guiPlugin::whatsThis() const
{
    return QString();
}

bool noise_marking_guiPlugin::isContainer() const
{
    return false;
}

QString noise_marking_guiPlugin::domXml() const
{
    return "<widget class=\"noise_marking_gui\" name=\"noise_marking_gui\">\n"
        " <property name=\"geometry\">\n"
        "  <rect>\n"
        "   <x>0</x>\n"
        "   <y>0</y>\n"
        "   <width>100</width>\n"
        "   <height>100</height>\n"
        "  </rect>\n"
        " </property>\n"
        "</widget>\n";
}

QString noise_marking_guiPlugin::includeFile() const
{
    return "noise_marking_gui.h";
}
