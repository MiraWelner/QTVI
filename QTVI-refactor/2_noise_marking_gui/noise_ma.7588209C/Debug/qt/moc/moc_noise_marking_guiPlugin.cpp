/****************************************************************************
** Meta object code from reading C++ file 'noise_marking_guiPlugin.h'
**
** Created by: The Qt Meta Object Compiler version 69 (Qt 6.10.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../../../noise_marking_guiPlugin.h"
#include <QtCore/qmetatype.h>
#include <QtCore/qplugin.h>

#include <QtCore/qtmochelpers.h>

#include <memory>


#include <QtCore/qxptype_traits.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'noise_marking_guiPlugin.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 69
#error "This file was generated using the moc from 6.10.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

#ifndef Q_CONSTINIT
#define Q_CONSTINIT
#endif

QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
QT_WARNING_DISABLE_GCC("-Wuseless-cast")
namespace {
struct qt_meta_tag_ZN23noise_marking_guiPluginE_t {};
} // unnamed namespace

template <> constexpr inline auto noise_marking_guiPlugin::qt_create_metaobjectdata<qt_meta_tag_ZN23noise_marking_guiPluginE_t>()
{
    namespace QMC = QtMocConstants;
    QtMocHelpers::StringRefStorage qt_stringData {
        "noise_marking_guiPlugin"
    };

    QtMocHelpers::UintData qt_methods {
    };
    QtMocHelpers::UintData qt_properties {
    };
    QtMocHelpers::UintData qt_enums {
    };
    return QtMocHelpers::metaObjectData<noise_marking_guiPlugin, qt_meta_tag_ZN23noise_marking_guiPluginE_t>(QMC::MetaObjectFlag{}, qt_stringData,
            qt_methods, qt_properties, qt_enums);
}
Q_CONSTINIT const QMetaObject noise_marking_guiPlugin::staticMetaObject = { {
    QMetaObject::SuperData::link<QObject::staticMetaObject>(),
    qt_staticMetaObjectStaticContent<qt_meta_tag_ZN23noise_marking_guiPluginE_t>.stringdata,
    qt_staticMetaObjectStaticContent<qt_meta_tag_ZN23noise_marking_guiPluginE_t>.data,
    qt_static_metacall,
    nullptr,
    qt_staticMetaObjectRelocatingContent<qt_meta_tag_ZN23noise_marking_guiPluginE_t>.metaTypes,
    nullptr
} };

void noise_marking_guiPlugin::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    auto *_t = static_cast<noise_marking_guiPlugin *>(_o);
    (void)_t;
    (void)_c;
    (void)_id;
    (void)_a;
}

const QMetaObject *noise_marking_guiPlugin::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *noise_marking_guiPlugin::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_staticMetaObjectStaticContent<qt_meta_tag_ZN23noise_marking_guiPluginE_t>.strings))
        return static_cast<void*>(this);
    if (!strcmp(_clname, "QDesignerCustomWidgetInterface"))
        return static_cast< QDesignerCustomWidgetInterface*>(this);
    if (!strcmp(_clname, "org.qt-project.QDesignerCustomWidgetInterface"))
        return static_cast< QDesignerCustomWidgetInterface*>(this);
    return QObject::qt_metacast(_clname);
}

int noise_marking_guiPlugin::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    return _id;
}

#ifdef QT_MOC_EXPORT_PLUGIN_V2
static constexpr unsigned char qt_pluginMetaDataV2_noise_marking_guiPlugin[] = {
    0xbf, 
    // "IID"
    0x02,  0x78,  0x30,  'o',  'r',  'g',  '.',  'q', 
    't',  '-',  'p',  'r',  'o',  'j',  'e',  'c', 
    't',  '.',  'Q',  't',  '.',  'Q',  'D',  'e', 
    's',  'i',  'g',  'n',  'e',  'r',  'C',  'u', 
    's',  't',  'o',  'm',  'W',  'i',  'd',  'g', 
    'e',  't',  'I',  'n',  't',  'e',  'r',  'f', 
    'a',  'c',  'e', 
    // "className"
    0x03,  0x77,  'n',  'o',  'i',  's',  'e',  '_', 
    'm',  'a',  'r',  'k',  'i',  'n',  'g',  '_', 
    'g',  'u',  'i',  'P',  'l',  'u',  'g',  'i', 
    'n', 
    0xff, 
};
QT_MOC_EXPORT_PLUGIN_V2(noise_marking_guiPlugin, noise_marking_guiPlugin, qt_pluginMetaDataV2_noise_marking_guiPlugin)
#else
QT_PLUGIN_METADATA_SECTION
Q_CONSTINIT static constexpr unsigned char qt_pluginMetaData_noise_marking_guiPlugin[] = {
    'Q', 'T', 'M', 'E', 'T', 'A', 'D', 'A', 'T', 'A', ' ', '!',
    // metadata version, Qt version, architectural requirements
    0, QT_VERSION_MAJOR, QT_VERSION_MINOR, qPluginArchRequirements(),
    0xbf, 
    // "IID"
    0x02,  0x78,  0x30,  'o',  'r',  'g',  '.',  'q', 
    't',  '-',  'p',  'r',  'o',  'j',  'e',  'c', 
    't',  '.',  'Q',  't',  '.',  'Q',  'D',  'e', 
    's',  'i',  'g',  'n',  'e',  'r',  'C',  'u', 
    's',  't',  'o',  'm',  'W',  'i',  'd',  'g', 
    'e',  't',  'I',  'n',  't',  'e',  'r',  'f', 
    'a',  'c',  'e', 
    // "className"
    0x03,  0x77,  'n',  'o',  'i',  's',  'e',  '_', 
    'm',  'a',  'r',  'k',  'i',  'n',  'g',  '_', 
    'g',  'u',  'i',  'P',  'l',  'u',  'g',  'i', 
    'n', 
    0xff, 
};
QT_MOC_EXPORT_PLUGIN(noise_marking_guiPlugin, noise_marking_guiPlugin)
#endif  // QT_MOC_EXPORT_PLUGIN_V2

QT_WARNING_POP
