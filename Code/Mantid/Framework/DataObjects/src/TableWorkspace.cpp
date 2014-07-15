#include "MantidDataObjects/TableWorkspace.h"
#include "MantidKernel/Logger.h"
#include "MantidAPI/ColumnFactory.h"
#include "MantidAPI/WorkspaceProperty.h"
#include "MantidAPI/WorkspaceFactory.h"

#include <iostream>

namespace Mantid
{
  namespace DataObjects
  {
    namespace
    {
      /// static logger
      Kernel::Logger g_log("TableWorkspace");
    }

    DECLARE_WORKSPACE(TableWorkspace)


    /// Constructor
    TableWorkspace::TableWorkspace(size_t nrows) : ITableWorkspace(), m_rowCount(0),
    m_LogManager(new API::LogManager)
    {
      setRowCount(nrows);
    }
 

    ///Destructor
    TableWorkspace::~TableWorkspace()
    {}

    size_t TableWorkspace::getMemorySize() const
    {
      size_t data_size = 0;
      for(column_const_it c = m_columns.begin();c!=m_columns.end();c++)
      {
        data_size += (*c)->sizeOfData();

      }
      data_size+= m_LogManager->getMemorySize();
      return data_size;
    }

    /** @param type :: Data type of the column.
        @param name :: Column name.
        @return A shared pointer to the created column (will be null on failure).
    */
    API::Column_sptr TableWorkspace::addColumn(const std::string& type, const std::string& name)
    {
        API::Column_sptr c;
        if (type.empty())
        {
            g_log.error("Empty string passed as type argument of createColumn.");
            return c;
        }
        if (name.empty())
        {
            g_log.error("Empty string passed as name argument of createColumn.");
            return c;
        }
        // Check that there is no column with the same name.
        column_it ci = std::find_if(m_columns.begin(),m_columns.end(),FindName(name));
        if (ci != m_columns.end())
        {
            g_log.error()<<"Column with name "<<name<<" already exists.\n";
            return c;
        }
        try
        {
            c = API::ColumnFactory::Instance().create(type);
            m_columns.push_back(c);
            c->setName(name);
            resizeColumn(c.get(),rowCount());
        }
        catch(Kernel::Exception::NotFoundError& e)
        {
            g_log.error()<<"Column of type "<<type<<" and name "<<name<<" has not been created.\n";
            g_log.error()<<e.what()<<'\n';
            return c;
        }
        return c;
    }

    /** If count is greater than the current number of rows extra rows are added to the bottom of the table.
        Otherwise rows at the end are erased to reach the new size.
        @param count :: New number of rows.
    */
    void TableWorkspace::setRowCount(size_t count)
    {
        if (count == rowCount()) return;
        for(column_it ci=m_columns.begin();ci!=m_columns.end();ci++)
            resizeColumn(ci->get(),count);
        m_rowCount = count;
    }

    /// Gets the shared pointer to a column.
    API::Column_sptr TableWorkspace::getColumn(const std::string& name)
    {
        column_it ci = std::find_if(m_columns.begin(),m_columns.end(),FindName(name));
        if (ci == m_columns.end())
        {
            std::string str = "Column " + name + " does not exist.\n";
            g_log.error(str);
            throw std::runtime_error(str);
        }
        return *ci;
    }
   /// Gets the shared pointer to a column.
    API::Column_const_sptr TableWorkspace::getColumn(const std::string& name)const
    {
        column_const_it c_it  = m_columns.begin();
        column_const_it c_end = m_columns.end();
        for(;c_it!=c_end;c_it++){
            if(c_it->get()->name() == name){
                 return *c_it;
            }
        }
        std::string str = "Column " + name + " does not exist.\n";
        g_log.error(str);
        throw std::runtime_error(str);
    }

    /// Gets the shared pointer to a column.
    API::Column_sptr TableWorkspace::getColumn(size_t index)
    {
        if (index >= columnCount())
        {
            std::string str = "Column index is out of range";
            g_log.error()<<str<<": "<<index<<"("<<columnCount()<<")\n";
            throw std::range_error(str);
        }
        return m_columns[index];
    }

    /// Gets the shared pointer to a column.
    API::Column_const_sptr TableWorkspace::getColumn(size_t index) const
    {
        if (index >= columnCount())
        {
            std::string str = "Column index is out of range";
            g_log.error()<<str<<": "<<index<<"("<<columnCount()<<")\n";
            throw std::range_error(str);
        }
        return m_columns[index];
    }

    void TableWorkspace::removeColumn( const std::string& name)
    {
        column_it ci = std::find_if(m_columns.begin(),m_columns.end(),FindName(name));
        if (ci != m_columns.end())
        {
            if ( !ci->unique() )
            {
                g_log.error()<<"Deleting column in use ("<<name<<").\n";
            }
            m_columns.erase(ci);
        }
    }

    /** @param index :: Points where to insert the new row.
        @return Position of the inserted row.
    */
    size_t TableWorkspace::insertRow(size_t index)
    {
        if (index >= rowCount()) index = rowCount();
        for(column_it ci=m_columns.begin();ci!=m_columns.end();ci++)
            insertInColumn(ci->get(),index);
        ++m_rowCount;
        return index;
    }

    /** @param index :: Row to delete.
    */
    void TableWorkspace::removeRow(size_t index)
    {
        if (index >= rowCount())
        {
            g_log.error()<<"Attempt to delete a non-existing row ("<<index<<")\n";
            return;
        }
        for(column_it ci=m_columns.begin();ci!=m_columns.end();ci++)
            removeFromColumn(ci->get(),index);
        --m_rowCount;
    }

    std::vector<std::string> TableWorkspace::getColumnNames() const
    {
        std::vector<std::string> nameList;
        nameList.reserve(m_columns.size());
        for(auto ci=m_columns.begin();ci!=m_columns.end();ci++)
            nameList.push_back((*ci)->name());
        return nameList;
    }

    bool TableWorkspace::addColumn(boost::shared_ptr<API::Column> column)
    {
      column_it ci = std::find_if(m_columns.begin(),m_columns.end(),FindName(column->name()));
      if (ci != m_columns.end())
      {
        g_log.error()<<"Column with name "<<column->name()<<" already exists.\n";
        return false;
      }
      else
      {
        m_columns.push_back(column);
      }
      return true;
    }

    TableWorkspace* TableWorkspace::clone() const
    {
      TableWorkspace* copy = new TableWorkspace(this->m_rowCount);
      column_const_it it=m_columns.begin();
      while(it != m_columns.end())
      {
        copy->addColumn(boost::shared_ptr<API::Column>((*it)->clone()));
        it++;
      }
      // copy logs/properties.
      copy->m_LogManager = boost::make_shared<API::LogManager>(*this->m_LogManager);
      return copy;
    };

//    template<>
//    boost::tuples::null_type TableWorkspace::make_TupleRef< boost::tuples::null_type >(size_t j,const std::vector<std::string>& names,size_t i)
//    {return boost::tuples::null_type();}

  } // namespace DataObjects
} // namespace Mantid

///\cond TEMPLATE
namespace Mantid
{
  namespace Kernel
  {
    template<> DLLExport
    DataObjects::TableWorkspace_sptr IPropertyManager::getValue<DataObjects::TableWorkspace_sptr>(const std::string &name) const
    {
      PropertyWithValue<DataObjects::TableWorkspace_sptr>* prop =
        dynamic_cast<PropertyWithValue<DataObjects::TableWorkspace_sptr>*>(getPointerToProperty(name));
      if (prop)
      {
        return *prop;
      }
      else
      {
        std::string message = "Attempt to assign property "+ name +" to incorrect type. Expected TableWorkspace.";
        throw std::runtime_error(message);
      }
    }

  } // namespace Kernel
} // namespace Mantid

///\endcond TEMPLATE
