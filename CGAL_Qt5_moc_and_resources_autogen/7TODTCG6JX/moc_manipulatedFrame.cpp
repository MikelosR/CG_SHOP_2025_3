/****************************************************************************
** Meta object code from reading C++ file 'manipulatedFrame.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.15.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include <memory>
#include "../../../../CGAL-5.6.1/include/CGAL/Qt/manipulatedFrame.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'manipulatedFrame.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.15.3. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_CGAL__qglviewer__ManipulatedFrame_t {
    QByteArrayData data[18];
    char stringdata0[285];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_CGAL__qglviewer__ManipulatedFrame_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_CGAL__qglviewer__ManipulatedFrame_t qt_meta_stringdata_CGAL__qglviewer__ManipulatedFrame = {
    {
QT_MOC_LITERAL(0, 0, 33), // "CGAL::qglviewer::ManipulatedF..."
QT_MOC_LITERAL(1, 34, 11), // "manipulated"
QT_MOC_LITERAL(2, 46, 0), // ""
QT_MOC_LITERAL(3, 47, 4), // "spun"
QT_MOC_LITERAL(4, 52, 22), // "setRotationSensitivity"
QT_MOC_LITERAL(5, 75, 11), // "sensitivity"
QT_MOC_LITERAL(6, 87, 25), // "setTranslationSensitivity"
QT_MOC_LITERAL(7, 113, 22), // "setSpinningSensitivity"
QT_MOC_LITERAL(8, 136, 19), // "setWheelSensitivity"
QT_MOC_LITERAL(9, 156, 18), // "setZoomSensitivity"
QT_MOC_LITERAL(10, 175, 21), // "setSpinningQuaternion"
QT_MOC_LITERAL(11, 197, 10), // "Quaternion"
QT_MOC_LITERAL(12, 208, 18), // "spinningQuaternion"
QT_MOC_LITERAL(13, 227, 13), // "startSpinning"
QT_MOC_LITERAL(14, 241, 14), // "updateInterval"
QT_MOC_LITERAL(15, 256, 12), // "stopSpinning"
QT_MOC_LITERAL(16, 269, 4), // "spin"
QT_MOC_LITERAL(17, 274, 10) // "spinUpdate"

    },
    "CGAL::qglviewer::ManipulatedFrame\0"
    "manipulated\0\0spun\0setRotationSensitivity\0"
    "sensitivity\0setTranslationSensitivity\0"
    "setSpinningSensitivity\0setWheelSensitivity\0"
    "setZoomSensitivity\0setSpinningQuaternion\0"
    "Quaternion\0spinningQuaternion\0"
    "startSpinning\0updateInterval\0stopSpinning\0"
    "spin\0spinUpdate"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_CGAL__qglviewer__ManipulatedFrame[] = {

 // content:
       8,       // revision
       0,       // classname
       0,    0, // classinfo
      12,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       2,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    0,   74,    2, 0x06 /* Public */,
       3,    0,   75,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
       4,    1,   76,    2, 0x0a /* Public */,
       6,    1,   79,    2, 0x0a /* Public */,
       7,    1,   82,    2, 0x0a /* Public */,
       8,    1,   85,    2, 0x0a /* Public */,
       9,    1,   88,    2, 0x0a /* Public */,
      10,    1,   91,    2, 0x0a /* Public */,
      13,    1,   94,    2, 0x0a /* Public */,
      15,    0,   97,    2, 0x0a /* Public */,
      16,    0,   98,    2, 0x09 /* Protected */,
      17,    0,   99,    2, 0x08 /* Private */,

 // signals: parameters
    QMetaType::Void,
    QMetaType::Void,

 // slots: parameters
    QMetaType::Void, QMetaType::QReal,    5,
    QMetaType::Void, QMetaType::QReal,    5,
    QMetaType::Void, QMetaType::QReal,    5,
    QMetaType::Void, QMetaType::QReal,    5,
    QMetaType::Void, QMetaType::QReal,    5,
    QMetaType::Void, 0x80000000 | 11,   12,
    QMetaType::Void, QMetaType::Int,   14,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void CGAL::qglviewer::ManipulatedFrame::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        auto *_t = static_cast<ManipulatedFrame *>(_o);
        (void)_t;
        switch (_id) {
        case 0: _t->manipulated(); break;
        case 1: _t->spun(); break;
        case 2: _t->setRotationSensitivity((*reinterpret_cast< qreal(*)>(_a[1]))); break;
        case 3: _t->setTranslationSensitivity((*reinterpret_cast< qreal(*)>(_a[1]))); break;
        case 4: _t->setSpinningSensitivity((*reinterpret_cast< qreal(*)>(_a[1]))); break;
        case 5: _t->setWheelSensitivity((*reinterpret_cast< qreal(*)>(_a[1]))); break;
        case 6: _t->setZoomSensitivity((*reinterpret_cast< qreal(*)>(_a[1]))); break;
        case 7: _t->setSpinningQuaternion((*reinterpret_cast< const Quaternion(*)>(_a[1]))); break;
        case 8: _t->startSpinning((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 9: _t->stopSpinning(); break;
        case 10: _t->spin(); break;
        case 11: _t->spinUpdate(); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        {
            using _t = void (ManipulatedFrame::*)();
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&ManipulatedFrame::manipulated)) {
                *result = 0;
                return;
            }
        }
        {
            using _t = void (ManipulatedFrame::*)();
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&ManipulatedFrame::spun)) {
                *result = 1;
                return;
            }
        }
    }
}

QT_INIT_METAOBJECT const QMetaObject CGAL::qglviewer::ManipulatedFrame::staticMetaObject = { {
    QMetaObject::SuperData::link<Frame::staticMetaObject>(),
    qt_meta_stringdata_CGAL__qglviewer__ManipulatedFrame.data,
    qt_meta_data_CGAL__qglviewer__ManipulatedFrame,
    qt_static_metacall,
    nullptr,
    nullptr
} };


const QMetaObject *CGAL::qglviewer::ManipulatedFrame::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *CGAL::qglviewer::ManipulatedFrame::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_CGAL__qglviewer__ManipulatedFrame.stringdata0))
        return static_cast<void*>(this);
    if (!strcmp(_clname, "MouseGrabber"))
        return static_cast< MouseGrabber*>(this);
    return Frame::qt_metacast(_clname);
}

int CGAL::qglviewer::ManipulatedFrame::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = Frame::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 12)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 12;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 12)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 12;
    }
    return _id;
}

// SIGNAL 0
void CGAL::qglviewer::ManipulatedFrame::manipulated()
{
    QMetaObject::activate(this, &staticMetaObject, 0, nullptr);
}

// SIGNAL 1
void CGAL::qglviewer::ManipulatedFrame::spun()
{
    QMetaObject::activate(this, &staticMetaObject, 1, nullptr);
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
