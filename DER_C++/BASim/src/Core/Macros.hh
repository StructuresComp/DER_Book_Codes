/**
 * \file Macros.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 08/31/2009
 */

#ifndef MACROS_HH
#define MACROS_HH

namespace BASim {

#define BA_ADD_PROPERTY(handle, container, size)                  \
  template <typename T> inline void                               \
  add_property(handle<T>& ph, const std::string& name, T t = T()) \
  {                                                               \
    ph = handle<T>( container.add<T>(name) );                     \
    container.resize(size);                                       \
    container.property(ph).set_default(t);                        \
  }

#define BA_ACCESS_PROPERTY(handle, container)     \
  template <typename T> inline Property<T>&       \
  property(const handle<T>& ph)                   \
  {                                               \
    return container.property(ph);                \
  }                                               \
                                                  \
  template <typename T> inline const Property<T>& \
  property(const handle<T>& ph) const             \
  {                                               \
    return container.property(ph);                \
  }

#define BA_ACCESS_SINGULAR_PROPERTY(handle, container)                \
  template <typename T> inline typename Property<T>::reference        \
  property(const handle<T>& ph)                                       \
  {                                                                   \
    return container.property(ph)[0];                                 \
  }                                                                   \
                                                                      \
  template <typename T> inline typename Property<T>::const_reference  \
  property(const handle<T>& ph) const                                 \
  {                                                                   \
    return container.property(ph)[0];                                 \
  }

#define BA_PROPERTY_EXISTS(handle, container)                       \
  template <typename T> inline bool                                 \
  property_exists(const handle<T>&, const std::string& name) const  \
  {                                                                 \
    return container.exists<T>(name);                               \
  }

#define BA_PROPERTY_HANDLE(hndl, container, size)                   \
  template <typename T> inline void                                 \
  property_handle(hndl<T>& ph, const std::string& name, T t = T())  \
  {                                                                 \
    if ( container.exists<T>(name) ) {                              \
      ph = hndl<T>( container.handle<T>(name) );                    \
    } else {                                                        \
      ph = hndl<T>( container.add<T>(name) );                       \
      container.resize(size);                                       \
      container.property(ph).set_default(t);                        \
    }                                                               \
  }

#define BA_ACCESS_CONTAINER(handle, cont)         \
  template <typename T> PropertyContainer&        \
  container(const handle<T>&)                     \
  {                                               \
    return cont;                                  \
  }                                               \
                                                  \
  template <typename T> const PropertyContainer&  \
  container(const handle<T>&) const               \
  {                                               \
    return cont;                                  \
  }

#define BA_CREATE_PROPERTY(handle, container, size) \
  BA_ADD_PROPERTY(handle, container, size)          \
  BA_ACCESS_PROPERTY(handle, container)             \
  BA_PROPERTY_EXISTS(handle, container)             \
  BA_PROPERTY_HANDLE(handle, container, size)       \
  BA_ACCESS_CONTAINER(handle, container)

#define BA_CREATE_SINGULAR_PROPERTY(handle, container)  \
  BA_ADD_PROPERTY(handle, container, 1)                 \
  BA_ACCESS_SINGULAR_PROPERTY(handle, container)        \
  BA_PROPERTY_EXISTS(handle, container)                 \
  BA_PROPERTY_HANDLE(handle, container, 1)              \
  BA_ACCESS_CONTAINER(handle, container)

#define BA_INHERIT_BASE(base)                   \
  using base::add_property;                     \
  using base::property;                         \
  using base::property_exists;                  \
  using base::property_handle;                  \
  using base::container;

} // namespace BASim

#endif // MACROS_HH
